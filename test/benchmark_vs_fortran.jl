# Kraken.jl against Fortran KRAKEN: forward runs and gradients
# =============================================================
#
# WHAT THIS MEASURES
#   1. The forward solve. `kraken_jl` against unmodified Fortran `kraken.exe` from
#      AcousticsToolbox_jll, on the same environments, after checking that the two agree.
#   2. Gradients of a scalar loss over the wavenumbers, five ways: forward mode (ForwardDiff),
#      reverse mode (Zygote and Mooncake), central finite differences of `kraken_jl`, and central
#      finite differences of `kraken.exe` itself (2N process launches).
#   3. How the gradient cost grows with the number of parameters N, using an N-point sound-speed
#      profile as the unknown.
#
# HOW TO RUN IT
#   It runs in the test environment, which already has every package used here plus the Fortran
#   harness in `test/reference/`. Once per clone or worktree, dev-link that environment so it loads
#   this checkout's Kraken rather than the registered release:
#
#       julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#
#   then, from the repository root:
#
#       julia --project=test test/benchmark_vs_fortran.jl
#
#   It is written in `#%%` cells, so it also runs block by block in VS Code. Run the cells in order;
#   each one only needs the cells above it.
#
#   Knobs, all environment variables:
#     KRAKEN_BENCH_SECONDS   time budget per timed quantity (default 2). Raise it for steadier medians.
#     KRAKEN_BENCH_OUT       a directory; if set, the last cell writes the three tables there as CSV.
#     KRAKEN_FORTRAN_BIN     a local Acoustics Toolbox `bin/` to benchmark instead of the jll binaries.
#
# READING THE NUMBERS FAIRLY
#   * `kraken.exe` is a separate process. Every Fortran call pays a process launch, writes a `.env`
#     and reads a `.mod`/`.prt` back. Section 3 reports the binary on its own, the floor of
#     launching any process at all, and the full harness round trip, so you can see how much of the
#     Fortran time is the solve.
#   * Both solvers pick their own mesh and Richardson-extrapolate across refinements, but they do not
#     pick the same meshes, so the forward times compare two complete solves of the same problem,
#     not two identical computations.
#   * Times are minimum and median over repeated runs after a warm-up call, so compilation is
#     excluded. Section 4 reports the first-call time separately; for Mooncake that is mostly
#     building its rules, and it is paid once per session.

#%% ---------------------------------------------------------------------------------------------
# 1. Setup
# -----------------------------------------------------------------------------------------------

using Kraken
using BenchmarkTools
using DifferentiationInterface
using Printf
using Statistics
using FiniteDiff: FiniteDiff
using ForwardDiff: ForwardDiff
using Mooncake: Mooncake
using Zygote: Zygote

include(joinpath(@__DIR__, "reference", "KrakenReference.jl"))
const KR = KrakenReference

KR.fortran_available() ||
    error("No Fortran KRAKEN binary found. Check AcousticsToolbox_jll, or set KRAKEN_FORTRAN_BIN.")

const BUDGET = parse(Float64, get(ENV, "KRAKEN_BENCH_SECONDS", "2"))

# Root tolerances for every gradient below. They equal `kraken_jl`'s default and are spelled out so
# the gradients do not move if it changes: at 1e-6 the coarse-mesh roots carry enough truncation to
# shift the derivative by a few 1e-6 (the TOL docstring in test/reverse_ad_tests.jl). The forward
# benchmark uses the defaults, since that is what a user runs.
const TOL = (abstol=1e-10, reltol=1e-10)

println("Kraken.jl  ", pkgversion(Kraken), "  at ", pathof(Kraken))
println("Fortran    ", KR.describe())
println("Julia      ", VERSION, ", ", Threads.nthreads(), " thread(s), ", Sys.cpu_info()[1].model)
println("Budget     ", BUDGET, " s per timed quantity")

"""
    bench(f; seconds=BUDGET) -> (min, med, n)

Time `f()` in seconds: one warm-up call, then repeated single evaluations for `seconds`. A call slower
than the budget still gets one timed sample.
"""
function bench(f; seconds=BUDGET)
    f()
    trial = run(@benchmarkable($f(), evals = 1); seconds=seconds, samples=100_000)
    return (min=minimum(trial.times) / 1e9, med=median(trial.times) / 1e9, n=length(trial.times))
end

# Relative difference, with a floor so an entry near zero is not judged against its own rounding.
relerr(a, b) = abs(a - b) / max(abs(b), 1e-12)

ms(t) = 1e3 * t

#%% ---------------------------------------------------------------------------------------------
# 2. Test cases, and a check that both solvers find the same modes
# -----------------------------------------------------------------------------------------------
# Timing two solvers is only meaningful if they solve the same problem, so this cell runs each case
# once on both sides and compares the wavenumbers. Fortran's are taken at the best precision it
# reports: the `.prt`'s ten printed digits where it lists the mode, the single-precision `.mod`
# otherwise. Agreement to ~1e-6 relative is what the cross-validation tests assert.

forward_cases = [
    (name="Pekeris 100 Hz", env=UnderwaterEnv(pekeris_env()...), freq=100.0),
    (name="Pekeris 500 Hz", env=UnderwaterEnv(pekeris_env()...), freq=500.0),
    (name="One layer 100 Hz", env=UnderwaterEnv(one_layer_env()...), freq=100.0),
    (name="Two sloped layers 100 Hz", env=UnderwaterEnv(two_layer_slope_env()...), freq=100.0),
    (name="Munk 50 Hz", env=UnderwaterEnv(munk_env()...), freq=50.0),
]

println("\n2. AGREEMENT (default kraken_jl tolerances, Fortran automatic mesh)")
@printf("   %-26s %8s %8s %14s\n", "case", "Julia M", "Fortran M", "max rel Δkr")
for c in forward_cases
    sol = kraken_jl(c.env, c.freq)
    kr_f = KR._best_fortran_kr(KR.run_fortran_kraken(c.env, c.freq))
    n = min(length(sol.kr), length(kr_f))
    worst = maximum(relerr(real(sol.kr[m]), kr_f[m]) for m in 1:n)
    @printf("   %-26s %8d %8d %14.2e\n", c.name, length(sol.kr), length(kr_f), worst)
end

#%% ---------------------------------------------------------------------------------------------
# 3. Forward run
# -----------------------------------------------------------------------------------------------
# Four times per case:
#   Julia         `kraken_jl(env, freq)`, in-process.
#   kraken.exe    the binary alone, on a `.env` written once beforehand.
#   spawn floor   launching `true`, a process that does nothing. Subtract it from the binary's time
#                 to estimate the solve alone.
#   harness       `run_fortran_kraken`: write the `.env`, launch, parse the `.mod` and `.prt`. This
#                 is what a Julia caller pays to get Fortran's answer.

"""
A zero-argument closure that runs `kraken.exe` on `env` written once to a fresh directory, and that
directory. Output goes to `/dev/null`, because the harness discards it too.
"""
function fortran_binary(env, freq; kwargs...)
    dir = mktempdir()
    KR.write_env_file(joinpath(dir, "case"), env, freq; kwargs...)
    cmd = Cmd(`$(KR.kraken_cmd()) case`; dir=dir)
    return () -> run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull)), dir
end

spawn_floor = bench(() -> run(`true`))

forward_results = map(forward_cases) do c
    t_julia = bench(() -> kraken_jl(c.env, c.freq))
    exe, dir = fortran_binary(c.env, c.freq)
    t_exe = bench(exe)
    rm(dir; recursive=true, force=true)
    t_harness = bench(() -> KR.run_fortran_kraken(c.env, c.freq))
    return (; c.name, nmodes=length(kraken_jl(c.env, c.freq).kr), t_julia, t_exe, t_harness)
end

println("\n3. FORWARD RUN (median ms; spawn floor $(@sprintf("%.2f", ms(spawn_floor.med))) ms)")
@printf(
    "   %-26s %5s %10s %10s %12s %10s %12s\n",
    "case",
    "M",
    "Julia",
    "kraken.exe",
    "exe − floor",
    "harness",
    "exe / Julia"
)
for r in forward_results
    @printf(
        "   %-26s %5d %10.3f %10.3f %12.3f %10.3f %12.2f\n",
        r.name,
        r.nmodes,
        ms(r.t_julia.med),
        ms(r.t_exe.med),
        ms(r.t_exe.med - spawn_floor.med),
        ms(r.t_harness.med),
        r.t_exe.med / r.t_julia.med
    )
end

#%% ---------------------------------------------------------------------------------------------
# 4. Gradients: forward mode, reverse mode, and finite differences on both sides
# -----------------------------------------------------------------------------------------------
# The loss is L(θ) = Σₘ kᵣₘ over the first K modes, through `kraken_jl` (mesh refinement and
# Richardson extrapolation included). K stops one short of the modes both solvers find, so a ±h
# perturbation cannot push the last counted mode past cutoff, and K ≤ 30 so every counted mode is
# one the `.prt` prints to ten digits.
#
# FINITE DIFFERENCES OF FORTRAN
#   Central differences, 2N launches of `kraken.exe` at a relative step `h = 1e-3`. The step is
#   large because the `.prt` prints ten digits of kᵣ ~ 0.4; at 1e-6 the difference would be mostly
#   rounding. Fortran's mesh is pinned to Julia's coarsest one (`props.Nz_vec` at θ₀), so both
#   perturbed runs solve on the same grid; KRAKEN's automatic mesh would move with cb and the layer
#   thicknesses and the difference would measure the move.
#
# WHICH PARAMETERS FINITE DIFFERENCES CAN MEASURE
#   Julia's mesh point count is `ceil` of something depending on cb and every layer thickness, so
#   `kraken_jl` jumps where it changes. Both AD modes report the derivative at a fixed mesh, which is
#   the meaningful quantity; a finite difference of `kraken_jl` across a jump measures the jump. The
#   accuracy column marks those parameters with `*` and leaves them out of the worst-case figure.
#   Fortran FD is immune because its mesh is pinned.

pekeris_θ(θ) = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)

function onelayer_θ(θ)
    return UnderwaterEnv(one_layer_env(; c0=θ[1], c1=θ[2], cb=θ[3], ρ0=θ[4], ρ1=θ[5], ρb=θ[6], h0=θ[7], h1=θ[8])...)
end

# Munk as a function of the five numbers that define it, θ = [c_ref, ε, z_axis, cb, ρb]; at the
# values below it reproduces `munk_env()`. Copied from test/reverse_ad_tests.jl. The knot depths are
# built by comprehension, not `range`, because Zygote has no adjoint for a Float64 `range`.
function munk_θ(θ; depth=5000.0, nsamples=51)
    c_ref, ϵ, z_axis, cb, ρb = θ
    ρ0 = 1000.0
    zs = [depth * (k - 1) / (nsamples - 1) for k in 1:nsamples]
    profile(z) = (ẑ=2 * (z - z_axis) / z_axis; c_ref * (1 + ϵ * (ẑ - 1 + exp(-ẑ))))
    ssp = reduce(vcat, [[z profile(z) 0.0 ρ0 0.0 0.0] for z in zs])
    return UnderwaterEnv(ssp, [0.0 0.0 depth], [0.0 343.0 0.0 0.00121 0.0 0.0; depth cb 0.0 ρb 0.0 0.0])
end

# `fd_ok`: parameters a finite difference of `kraken_jl` measures (no mesh jump), from
# test/reverse_ad_tests.jl's ENVIRONMENTS table.
gradient_cases = [
    (
        name="Pekeris 100 Hz",
        envf=pekeris_θ,
        θ=[1500.0, 1600.0, 1000.0, 1500.0, 100.0],
        names=["c0", "cb", "ρ0", "ρb", "depth"],
        freq=100.0,
        fd_ok=[1, 3, 4],
    ),
    (
        name="One layer 100 Hz",
        envf=onelayer_θ,
        θ=[1500.0, 1550.0, 1600.0, 1000.0, 1500.0, 2000.0, 100.0, 20.0],
        names=["c0", "c1", "cb", "ρ0", "ρ1", "ρb", "h0", "h1"],
        freq=100.0,
        fd_ok=[1, 2, 4, 5, 6],
    ),
    (
        name="Munk 10 Hz",
        envf=munk_θ,
        θ=[1500.0, 0.00737, 1300.0, 1600.0, 1500.0],
        names=["c_ref", "ε", "z_axis", "cb", "ρb"],
        freq=10.0,
        fd_ok=[1, 2, 3, 5],
    ),
]

const BACKENDS = [
    ("ForwardDiff", AutoForwardDiff()),
    ("Zygote", AutoZygote()),
    ("Mooncake", AutoMooncake()),
    ("FiniteDiff (Julia)", AutoFiniteDiff(; fdtype=Val(:central))),
]

"""
    gradient_setup(case) -> (; loss, fortran_loss, K, nmesh)

The Julia loss, its Fortran counterpart on a pinned mesh, the mode count K they sum over, and that
mesh.
"""
function gradient_setup(case)
    env0 = case.envf(case.θ)
    sol0 = kraken_jl(env0, case.freq; TOL...)
    nmesh = AcousticProblemProperties(env0, case.freq).Nz_vec
    nf = length(KR._best_fortran_kr(KR.run_fortran_kraken(env0, case.freq; nmesh=nmesh)))
    K = clamp(min(length(sol0.kr), nf) - 1, 1, 30)
    loss(θ) = sum(kraken_jl(case.envf(θ), case.freq; TOL...).kr[1:K])
    fortran_loss(θ) = sum(KR._best_fortran_kr(KR.run_fortran_kraken(case.envf(θ), case.freq; nmesh=nmesh))[1:K])
    return (; loss, fortran_loss, K, nmesh)
end

"""Central differences of `f` at relative step `h`: 2N calls."""
function fd_relative(f, θ; h=1e-3)
    return map(eachindex(θ)) do k
        θp = copy(θ)
        θm = copy(θ)
        θp[k] *= 1 + h
        θm[k] *= 1 - h
        return (f(θp) - f(θm)) / (2 * θ[k] * h)
    end
end

gradient_results = []
for case in gradient_cases
    s = gradient_setup(case)
    N = length(case.θ)
    rows = []

    t_primal = bench(() -> s.loss(case.θ))
    push!(rows, (method="primal (Julia)", g=nothing, first=NaN, t=t_primal))

    for (label, backend) in BACKENDS
        local prep, g
        t_first = @elapsed begin
            prep = prepare_gradient(s.loss, backend, case.θ)
            g = gradient(s.loss, prep, backend, case.θ)
        end
        t = bench(() -> gradient(s.loss, prep, backend, case.θ))
        push!(rows, (method=label, g=g, first=t_first, t=t))
    end

    t_fprimal = bench(() -> s.fortran_loss(case.θ))
    push!(rows, (method="primal (Fortran)", g=nothing, first=NaN, t=t_fprimal))
    g_ffd = fd_relative(s.fortran_loss, case.θ)
    t_ffd = bench(() -> fd_relative(s.fortran_loss, case.θ); seconds=min(BUDGET, 1.0))
    push!(rows, (method="FiniteDiff (Fortran)", g=g_ffd, first=NaN, t=t_ffd))

    push!(gradient_results, (; case.name, N, s.K, s.nmesh, rows))

    # Reference for the accuracy column: Zygote, which ForwardDiff and Mooncake must match to
    # ~1e-9. The two finite differences are judged against it on the parameters they can measure.
    g_ref = only(r.g for r in rows if r.method == "Zygote")

    println("\n4. GRADIENT, $(case.name): N = $N parameters, L = Σ kᵣ over K = $(s.K) modes, mesh $(s.nmesh)")
    @printf("   %-22s %12s %12s %10s %14s\n", "method", "median ms", "× primal", "first s", "max rel err")
    for r in rows
        err = if r.g === nothing
            ""
        else
            idx = r.method == "FiniteDiff (Julia)" ? case.fd_ok : eachindex(case.θ)
            @sprintf("%.1e", maximum(relerr(r.g[i], g_ref[i]) for i in idx))
        end
        base = startswith(r.method, "FiniteDiff (Fortran)") || r.method == "primal (Fortran)" ? t_fprimal : t_primal
        @printf(
            "   %-22s %12.3f %12.1f %10s %14s\n",
            r.method,
            ms(r.t.med),
            r.t.med / base.med,
            isnan(r.first) ? "" : @sprintf("%.2f", r.first),
            err
        )
    end
    println("   (× primal is against the Julia primal, or the Fortran one for the Fortran rows)")

    # The gradient itself, so a large error can be traced to the parameter responsible.
    @printf("\n   %-8s %14s %14s %14s %14s\n", "param", "Zygote", "ForwardDiff", "FD (Julia)", "FD (Fortran)")
    pick(m) = only(r.g for r in rows if r.method == m)
    for i in eachindex(case.θ)
        mark = i in case.fd_ok ? " " : "*"
        @printf(
            "   %-8s %14.6e %14.6e %13.6e%s %14.6e\n",
            case.names[i],
            g_ref[i],
            pick("ForwardDiff")[i],
            pick("FiniteDiff (Julia)")[i],
            mark,
            g_ffd[i]
        )
    end
    println("   * mesh point count changes inside the Julia FD step; that entry measures the jump.")
end

#%% ---------------------------------------------------------------------------------------------
# 5. Cost against the number of parameters
# -----------------------------------------------------------------------------------------------
# The unknown is an N-point sound-speed profile in 100 m of water over a 1600 m/s half-space, at
# 100 Hz. Forward mode and both finite differences need work proportional to N; reverse mode needs
# one forward and one backward pass whatever N is. The depth knots are fixed outside the
# differentiated function.

function profile_env(cvec, z; cb=1600.0, ρ0=1000.0, ρb=1500.0)
    N = length(cvec)
    depth = z[end]
    ssp = hcat(z, cvec, zero(cvec), fill(ρ0, N), zero(cvec), zero(cvec))
    sspHS = [
        0.0 343.0 0.0 0.00121 0.0 0.0
        depth cb 0.0 ρb 0.0 0.0
    ]
    return UnderwaterEnv(ssp, [0.0 0.0 depth], sspHS)
end

const SCALING_N = (5, 10, 25, 50, 100)
const SCALING_FREQ = 100.0

scaling_results = map(SCALING_N) do N
    z = [100.0 * (k - 1) / (N - 1) for k in 1:N]
    c = 1500.0 .+ 5 .* sin.(3 .* (0:(N - 1)) ./ (N - 1))
    s = gradient_setup((envf=θ -> profile_env(θ, z), θ=c, freq=SCALING_FREQ))
    times = Dict{String,Float64}("primal (Julia)" => bench(() -> s.loss(c)).med)
    for (label, backend) in BACKENDS
        prep = prepare_gradient(s.loss, backend, c)
        times[label] = bench(() -> gradient(s.loss, prep, backend, c)).med
    end
    times["primal (Fortran)"] = bench(() -> s.fortran_loss(c)).med
    times["FiniteDiff (Fortran)"] = bench(() -> fd_relative(s.fortran_loss, c); seconds=min(BUDGET, 1.0)).med
    return (; N, times)
end

scaling_methods = [
    "primal (Julia)",
    "ForwardDiff",
    "Zygote",
    "Mooncake",
    "FiniteDiff (Julia)",
    "primal (Fortran)",
    "FiniteDiff (Fortran)",
]

println("\n5. GRADIENT COST vs N (median ms), N-point profile at $(SCALING_FREQ) Hz")
@printf("   %-22s", "method")
foreach(r -> @printf(" %10s", "N = $(r.N)"), scaling_results)
println()
for m in scaling_methods
    @printf("   %-22s", m)
    foreach(r -> @printf(" %10.2f", ms(r.times[m])), scaling_results)
    println()
end
println("\n   Same table as a multiple of the Julia primal:")
for m in scaling_methods[2:end]
    @printf("   %-22s", m)
    foreach(r -> @printf(" %10.1f", r.times[m] / r.times["primal (Julia)"]), scaling_results)
    println()
end

#%% ---------------------------------------------------------------------------------------------
# 6. Optional: write the tables as CSV
# -----------------------------------------------------------------------------------------------
# Set KRAKEN_BENCH_OUT to a directory to keep the numbers, e.g. for plotting. Times are seconds.

if haskey(ENV, "KRAKEN_BENCH_OUT")
    out = mkpath(ENV["KRAKEN_BENCH_OUT"])
    open(joinpath(out, "forward.csv"), "w") do io
        println(io, "case,nmodes,julia_med,julia_min,exe_med,exe_min,harness_med,spawn_floor_med")
        for r in forward_results
            println(
                io,
                join(
                    (
                        r.name,
                        r.nmodes,
                        r.t_julia.med,
                        r.t_julia.min,
                        r.t_exe.med,
                        r.t_exe.min,
                        r.t_harness.med,
                        spawn_floor.med,
                    ),
                    ",",
                ),
            )
        end
    end
    open(joinpath(out, "gradients.csv"), "w") do io
        println(io, "case,N,K,method,med,min,first_call")
        for g in gradient_results, r in g.rows
            println(io, join((g.name, g.N, g.K, r.method, r.t.med, r.t.min, r.first), ","))
        end
    end
    open(joinpath(out, "scaling.csv"), "w") do io
        println(io, "N,method,med")
        for r in scaling_results, m in scaling_methods
            println(io, join((r.N, m, r.times[m]), ","))
        end
    end
    println("\nWrote forward.csv, gradients.csv and scaling.csv to ", out)
end
