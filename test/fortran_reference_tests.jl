using Test
using Kraken

include("reference/KrakenReference.jl")
const KR = KrakenReference

# Cross-validation against unmodified Fortran KRAKEN. Plan tasks 3.3-3.7 extend this file with the
# .mod/.prt readers, the runner, and the per-environment comparison suite; for now it covers binary
# resolution and the .env writer.
#
# Everything here is gated on KR.fortran_available(): a platform with no AcousticsToolbox_jll build
# and no KRAKEN_FORTRAN_BIN must skip with a message, never error.

@testset "Fortran reference" begin
    @testset "Binary resolution" begin
        # ENV is global, and the whole point of these tests is to toggle the override, so restore
        # whatever the caller had regardless of how we exit.
        saved = get(ENV, KR.BIN_ENV_VAR, nothing)
        restore!() = saved === nothing ? delete!(ENV, KR.BIN_ENV_VAR) : (ENV[KR.BIN_ENV_VAR] = saved)

        try
            delete!(ENV, KR.BIN_ENV_VAR)

            if !KR.fortran_available()
                @info """
                No Fortran KRAKEN available — skipping cross-validation.
                AcousticsToolbox_jll has no build for this platform. Set $(KR.BIN_ENV_VAR) to a
                directory containing kraken.exe to use a local Acoustics Toolbox build.
                """
                @test_skip KR.fortran_available()
            else
                @info "Fortran reference resolved -- $(KR.describe())"

                # Without an override: the jll.
                @test KR.binary_source() === :jll
                @test KR.binary_source(; complex=true) === :jll
                @test isfile(KR.binary_path())
                @test isfile(KR.binary_path(; complex=true))
                @test basename(KR.binary_path()) == "kraken.exe"
                @test basename(KR.binary_path(; complex=true)) == "krakenc.exe"
                @test KR.override_dir() === nothing

                # The resolved binary must actually run. A path that exists but cannot execute
                # would make fortran_available() a false positive and turn every downstream test
                # into a confusing failure.
                #
                # KRAKEN exits 0 even on a fatal error -- it prints "STOP Fatal Error: Check the
                # print file for details" to stderr and stops. That is why the runner in task 3.4
                # scans the .prt file rather than trusting the exit code.
                mktempdir() do dir
                    err = IOBuffer()
                    cmd = Cmd(`$(KR.kraken_cmd()) no_such_env_file`; dir=dir)
                    run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=err))
                    @test occursin("Fatal Error", String(take!(err)))
                end

                # With an override pointing at a real local build: the override wins.
                local_build = "/Users/arielv/programs/AcousticsToolboxOALIB/bin"
                if isfile(joinpath(local_build, "kraken.exe"))
                    ENV[KR.BIN_ENV_VAR] = local_build
                    @test KR.binary_source() === :override
                    @test KR.binary_path() == joinpath(local_build, "kraken.exe")
                    @test KR.binary_path(; complex=true) == joinpath(local_build, "krakenc.exe")
                    @test KR.fortran_available()
                else
                    @info "No local Acoustics Toolbox build at $local_build — override test skipped."
                end

                # A set-but-unusable override falls back to the jll instead of failing. The variable
                # exists to *prefer* a local build, not to make the suite fragile on a machine where
                # someone left it pointing at a stale path.
                ENV[KR.BIN_ENV_VAR] = joinpath(tempdir(), "definitely-not-a-toolbox-build")
                @test KR.binary_source() === :jll
                @test KR.fortran_available()

                # An empty value is treated as unset, not as "the current directory".
                ENV[KR.BIN_ENV_VAR] = ""
                @test KR.override_dir() === nothing
                @test KR.binary_source() === :jll
            end
        finally
            restore!()
        end
    end

    # --- .env writer -------------------------------------------------------------------------

    # The environments the writer has to cover, paired with the checked-in file that documents the
    # format for that shape. `nothing` means there is no companion file — those cases exist to
    # exercise a many-point profile and the broadband record block.
    env_writer_cases = [
        (name="pekeris", env=UnderwaterEnv(pekeris_env()...), freq=100.0, ref="Pekeris_AV.env"),
        (name="onelayer", env=UnderwaterEnv(one_layer_env()...), freq=100.0, ref="onelayer_AV.env"),
        (name="onelayer_slope", env=UnderwaterEnv(one_layer_slope_env()...), freq=100.0, ref="onelayer_slope_AV.env"),
        (name="twolayer_slope", env=UnderwaterEnv(two_layer_slope_env()...), freq=100.0, ref="twolayer_slope_AV.env"),
        (name="munk", env=UnderwaterEnv(munk_env()...), freq=50.0, ref=nothing),
        (name="pekeris_broadband", env=UnderwaterEnv(pekeris_env()...), freq=[50.0, 100.0, 200.0], ref=nothing),
    ]

    # Numeric tokens on one .env record: strip the `!` comment, the quoted strings (title,
    # top/bottom options) and the `/` terminator, then parse what is left.
    function env_numbers(line)
        text = replace(first(split(line, '!')), r"'[^']*'" => " ", '/' => ' ')
        return Float64[v for v in (tryparse(Float64, t) for t in split(text)) if v !== nothing]
    end

    # Split an .env file into the three groups the writer is responsible for getting right, without
    # reimplementing the full reader that task 3.7 will write. The bottom-option record is the
    # anchor: it is the first quoted record after the top options on line 4, and everything before
    # it is the environment proper (freq, nmedia, and the per-medium mesh + SSP rows).
    function env_sections(path)
        lines = filter(l -> !isempty(strip(first(split(l, '!')))), collect(eachline(path)))
        bot = findfirst(i -> i > 4 && occursin('\'', lines[i]), eachindex(lines))
        bot === nothing && error("No bottom-option record found in $path")
        return (
            environment=reduce(vcat, env_numbers.(lines[1:(bot - 1)]); init=Float64[]),
            halfspace=env_numbers(lines[bot + 1]),
            speed_limits=env_numbers(lines[bot + 2]),
        )
    end

    @testset "env writer" begin
        @testset "$(case.name)" for case in env_writer_cases
            text = KR.env_file_string(case.env, case.freq)

            # Structural invariants that hold for every environment, checked on the text itself so
            # a failure points at the record rather than at a Fortran error message.
            @test occursin("! NMEDIA", text)
            @test occursin("! CLOW  CHIGH (m/s)", text)
            @test occursin(r"'CVW", text)                          # C-linear SSP, vacuum surface
            @test count(l -> occursin("! NMESH", l), split(text, '\n')) == length(case.env.layer_depth)
            @test occursin("! NFREQ", text) == (case.freq isa AbstractVector)

            if case.ref !== nothing
                # The generated file and the checked-in one describe the same environment, so every
                # record that *defines* the environment has to agree numerically. The three places
                # they legitimately differ are all search/output settings rather than physics:
                # CLOW (1400 by hand vs 0.99*min(c) here), and the source/receiver depth vectors.
                mine = mktempdir(dir -> env_sections(KR.write_env_file(joinpath(dir, case.name), case.env, case.freq)))
                theirs = env_sections(joinpath(@__DIR__, "standard_envs", case.ref))

                @test mine.environment == theirs.environment
                # The checked-in halfspace rows stop after 5 values and let Fortran default ASB to
                # zero; this writer always emits all 6.
                @test mine.halfspace[1:5] == theirs.halfspace[1:5]
                @test mine.halfspace[6] == 0.0
                # CHIGH is the half-space sound speed in both; CLOW differs by design but must
                # still sit below every sound speed in the environment.
                @test mine.speed_limits[2] == theirs.speed_limits[2]
                @test mine.speed_limits[1] < minimum(case.env.c.c)
            end
        end
    end

    # Run `kraken.exe` (or `krakenc.exe`) on a checked-in .env and hand back the file root. Task
    # 3.4 replaces this with the real runner; here it only needs to be enough to produce a .mod and
    # a .prt to read. Note the .mod/.prt files were untracked in task 2.4 and are gitignored, so
    # they are regenerated rather than read from the tree.
    function run_checked_in_env(dir, name; complex=false, as=name)
        cp(joinpath(@__DIR__, "standard_envs", name * ".env"), joinpath(dir, as * ".env"); force=true)
        cmd = Cmd(`$(KR.kraken_cmd(; complex=complex)) $(as)`; dir=dir)
        run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull))
        return joinpath(dir, as)
    end

    if KR.fortran_available()
        @testset ".mod and .prt readers" begin
            @testset "Pekeris_AV" begin
                mktempdir() do dir
                    root = run_checked_in_env(dir, "Pekeris_AV")
                    modes = KR.read_mod_file(root * ".mod")
                    grp = KR.read_grp(root * ".prt")

                    # The two files describe the same solve, so they must agree on how many modes
                    # there are and on what their wavenumbers are.
                    @test modes.nmodes == 5
                    @test length(grp.m) == modes.nmodes
                    @test grp.m == 1:modes.nmodes
                    @test real(modes.kᵣ[1]) ≈ 0.4179 atol = 1e-4
                    @test modes.freq == 100.0
                    @test modes.freqs == [100.0]

                    # The .mod stores wavenumbers in single precision while the .prt prints ten
                    # digits of the double-precision value, so they agree only to ~1e-7.
                    @test maximum(abs.(real.(modes.kᵣ) .- real.(grp.kᵣ)) ./ abs.(real.(grp.kᵣ))) < 1e-6
                    @test all(iszero, imag.(modes.kᵣ))   # no attenuation in this environment
                    @test all(iszero, imag.(grp.kᵣ))

                    # `depths` is zTab -- the union of the .env source and receiver depths, which
                    # for this file is 100 receivers over 0..100 m plus one source at 25 m.
                    @test length(modes.depths) == size(modes.ϕ, 1)
                    @test modes.depths[1] == 0.0
                    @test modes.depths[end] ≈ 100.0 atol = 1e-3
                    @test issorted(modes.depths)

                    # Pressure-release surface: every mode vanishes at z = 0.
                    @test maximum(abs.(modes.ϕ[1, :])) < 1e-10
                    @test all(>(0), maximum(abs.(modes.ϕ); dims=1))

                    # Phase speeds are ω/kᵣ and must sit between the water and half-space speeds.
                    @test grp.phase_speed ≈ 2π * 100.0 ./ real.(grp.kᵣ) rtol = 1e-6
                    @test all(1500.0 .< grp.phase_speed .< 1600.0)
                end
            end

            @testset "broadband record stepping" begin
                mktempdir() do dir
                    root = run_checked_in_env(dir, "Pekeris_AV_BroadBand")
                    blocks = KR.read_grp_blocks(root * ".prt")
                    @test [b.freq for b in blocks] == [50.0, 100.0, 200.0, 300.0, 500.0]

                    counts = Int[]
                    for b in blocks
                        # Stepping to frequency i requires reading M for every block before it, so
                        # a wrong stride shows up as a garbage mode count or a read past EOF.
                        modes = KR.read_mod_file(root * ".mod"; freq=b.freq)
                        @test modes.freq == b.freq
                        @test modes.freqs == [50.0, 100.0, 200.0, 300.0, 500.0]
                        @test modes.nmodes == length(b.m)
                        @test maximum(abs.(real.(modes.kᵣ) .- real.(b.kᵣ)) ./ abs.(real.(b.kᵣ))) < 1e-6
                        # The higher-frequency blocks carry a small non-zero alpha, and the .prt
                        # prints it with G10.2 -- two significant digits. That, not the .mod, is
                        # the limiting precision here; see MOD_WAVENUMBER_DIGITS.
                        @test all(abs.(imag.(modes.kᵣ) .- imag.(b.kᵣ)) .<= 0.02 .* abs.(imag.(b.kᵣ)) .+ 1e-12)
                        @test maximum(abs.(modes.ϕ[1, :])) < 1e-10
                        push!(counts, modes.nmodes)
                    end
                    # More modes fit in the waveguide as the frequency rises.
                    @test issorted(counts)
                    @test counts == [2, 5, 10, 15, 24]
                end
            end

            @testset "many modes: subsampled .prt table" begin
                mktempdir() do dir
                    root = joinpath(dir, "munk")
                    KR.write_env_file(root, UnderwaterEnv(munk_env()...), 50.0)
                    run(pipeline(ignorestatus(Cmd(`$(KR.kraken_cmd()) munk`; dir=dir)); stdout=devnull, stderr=devnull))
                    modes = KR.read_mod_file(root * ".mod")
                    grp = KR.read_grp(root * ".prt")

                    # `DO mode = 1, M, MAX(1, M/30)` -- with ~100 modes the table lists every third
                    # one, so the row count is far below the mode count and the indices matter.
                    @test modes.nmodes > 30
                    @test length(grp.m) < modes.nmodes
                    @test grp.m[1] == 1
                    @test allunique(grp.m)
                    @test issorted(grp.m)
                    @test last(grp.m) <= modes.nmodes
                    @test allequal(diff(grp.m))
                    @test first(diff(grp.m)) == max(1, modes.nmodes ÷ 30)
                    # Indexing the .mod wavenumbers by the .prt's mode numbers must line them up.
                    @test maximum(abs.(real.(modes.kᵣ[grp.m]) .- real.(grp.kᵣ)) ./ abs.(real.(grp.kᵣ))) < 1e-6
                end
            end

            @testset "group speeds" begin
                mktempdir() do dir
                    # AcousticsToolbox_jll v2025.9's kraken.exe zeroes the Group Speed column while
                    # getting everything else right; its krakenc.exe does not. Anything downstream
                    # that needs group speeds has to know which binary it is talking to.
                    plain = KR.read_grp(run_checked_in_env(dir, "Pekeris_AV") * ".prt")
                    cplx = KR.read_grp(run_checked_in_env(dir, "Pekeris_AV"; complex=true, as="pek_c") * ".prt")

                    @test KR.has_group_speeds(cplx)
                    @test all(1400.0 .< cplx.v .< 1600.0)
                    # Group speed is below phase speed in a waveguide with a positive-gradient bottom.
                    @test all(cplx.v .< cplx.phase_speed)
                    # Same solve either way: the wavenumbers agree even though VG does not.
                    @test plain.kᵣ ≈ cplx.kᵣ rtol = 1e-9
                    if !KR.has_group_speeds(plain)
                        @info "kraken.exe ($(KR.binary_source())) reports no group speeds; " *
                            "krakenc.exe does. Group-speed comparisons must use complex=true."
                    end
                end
            end
        end

        @testset "runner" begin
            pekeris = UnderwaterEnv(pekeris_env()...)

            @testset "solves the canonical Pekeris case" begin
                result = KR.run_fortran_kraken(pekeris, 100.0)
                @test result.nmodes == 5
                @test real(result.kᵣ[1]) ≈ 0.4179 atol = 1e-4
                @test result.freq == 100.0
                @test result.dir === nothing            # cleaned up when keep_files is false
                @test isempty(result.warnings)
                @test length(result.grp.m) == result.nmodes
                @test real(result.grp.kᵣ[1]) ≈ 0.4179 atol = 1e-4
                @test size(result.ϕ, 2) == result.nmodes
                @test size(result.ϕ, 1) == length(result.depths)
            end

            @testset "keyword arguments reach the .env writer" begin
                # `rd` is the mode-shape grid, so setting it must change what comes back.
                result = KR.run_fortran_kraken(pekeris, 100.0; rd=range(0.0, 100.0; length=51))
                @test length(result.depths) == 52       # 51 receivers plus the source depth
                @test result.nmodes == 5
            end

            @testset "broadband returns one result per frequency" begin
                results = KR.run_fortran_kraken(pekeris, [50.0, 100.0, 200.0])
                @test length(results) == 3
                @test [r.freq for r in results] == [50.0, 100.0, 200.0]
                @test [r.nmodes for r in results] == [2, 5, 9]
                @test issorted([r.nmodes for r in results])
            end

            @testset "krakenc.exe via complex=true" begin
                result = KR.run_fortran_kraken(pekeris, 100.0; complex=true)
                @test result.nmodes == 5
                @test occursin("krakenc", result.binary)
                # This is the binary that actually reports group speeds -- see has_group_speeds.
                @test KR.has_group_speeds(result.grp)
                @test all(1400.0 .< result.grp.v .< 1600.0)
            end

            @testset "a failed run raises, quoting the Fortran" begin
                # NMESH=1 is below half the mesh KRAKEN wants, which ReadEnvironment rejects. The
                # binary still exits 0, so this is exactly the case an exit-code check would miss
                # and then misreport as a missing .mod file.
                err = try
                    KR.run_fortran_kraken(pekeris, 100.0; nmesh=1)
                    nothing
                catch e
                    e
                end
                @test err isa KR.FortranKrakenError
                message = sprint(showerror, err)
                @test occursin("FATAL ERROR", message)
                @test occursin("Mesh is too coarse", message)
                @test occursin("ReadEnvironment", message)
                @test occursin("kraken.exe", message)
                @test occursin("Mesh is too coarse", err.report)
            end

            # Run directories are `mktempdir()` children holding a `case.env`. Identifying them by
            # that file rather than by counting entries keeps these assertions honest when another
            # process is also writing to the system temp directory.
            run_dirs() = filter(readdir(tempdir(); join=true)) do path
                isdir(path) && isfile(joinpath(path, "case.env"))
            end

            @testset "temporary directories are cleaned up" begin
                before = Set(run_dirs())
                KR.run_fortran_kraken(pekeris, 100.0)
                for _ in 1:3
                    try
                        KR.run_fortran_kraken(pekeris, 100.0; nmesh=1)
                    catch
                    end
                end
                # Successful and failed runs alike must leave nothing behind.
                @test isempty(setdiff(Set(run_dirs()), before))
            end

            @testset "keep_files preserves the run directory" begin
                result = KR.run_fortran_kraken(pekeris, 100.0; keep_files=true)
                try
                    @test result.dir !== nothing
                    @test isdir(result.dir)
                    @test sort(readdir(result.dir)) == ["case.env", "case.mod", "case.prt"]
                finally
                    rm(result.dir; recursive=true, force=true)
                end

                # keep_files survives a failure too -- that is exactly when you want to look at the
                # inputs. The directory is not reachable through the exception, so find it by diff.
                before = Set(run_dirs())
                try
                    KR.run_fortran_kraken(pekeris, 100.0; nmesh=1, keep_files=true)
                catch
                end
                kept = collect(setdiff(Set(run_dirs()), before))
                @test length(kept) == 1
                for dir in kept
                    @test isfile(joinpath(dir, "case.env"))     # the input that caused the failure
                    @test isfile(joinpath(dir, "case.prt"))     # and the Fortran's own diagnosis
                    rm(dir; recursive=true, force=true)
                end
            end

            @testset "bindir overrides the binary" begin
                local_build = "/Users/arielv/programs/AcousticsToolboxOALIB/bin"
                if isfile(joinpath(local_build, "kraken.exe"))
                    result = KR.run_fortran_kraken(pekeris, 100.0; bindir=local_build)
                    @test result.binary == joinpath(local_build, "kraken.exe")
                    @test result.nmodes == 5
                    # The 2023 build does report group speeds where the jll's kraken.exe does not.
                    @test KR.has_group_speeds(result.grp)
                else
                    @info "No local Acoustics Toolbox build at $local_build — bindir test skipped."
                end
                @test_throws ErrorException KR.run_fortran_kraken(
                    pekeris, 100.0; bindir=joinpath(tempdir(), "not-a-toolbox")
                )
            end
        end

        @testset "env writer accepted by kraken.exe" begin
            @testset "$(case.name)" for case in env_writer_cases
                mktempdir() do dir
                    root = joinpath(dir, case.name)
                    KR.write_env_file(root, case.env, case.freq)
                    cmd = Cmd(`$(KR.kraken_cmd()) $(case.name)`; dir=dir)
                    run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull))

                    # kraken.exe exits 0 even on a fatal error, so the .prt is the only honest
                    # signal -- see the architecture note in the plan.
                    prt = root * ".prt"
                    @test isfile(prt)
                    report = read(prt, String)
                    bad = filter(l -> occursin("ERROR", uppercase(l)), split(report, '\n'))
                    isempty(bad) || @info "kraken.exe reported errors for $(case.name)" bad
                    @test isempty(bad)

                    # A .mod file that exists and is non-empty means it got all the way through the
                    # solve, not just through the reader.
                    @test isfile(root * ".mod")
                    @test filesize(root * ".mod") > 0
                    @test occursin("Number of modes", report)
                end
            end
        end
    end
end
