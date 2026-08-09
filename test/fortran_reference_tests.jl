using Test
using Kraken

# Task 4.7 validates the Milestone 4 gradients against kraken.exe, so this file needs both AD
# backends. Both are already test-environment dependencies for test/reverse_ad_tests.jl.
using ForwardDiff
using Zygote

include("reference/KrakenReference.jl")
const KR = KrakenReference

# Cross-validation against unmodified Fortran KRAKEN: binary resolution, the .env writer, the
# .mod/.prt readers, the runner, the comparison utility, and the per-environment regression sweep.
# Plan task 3.7 extends the sweep with cases parsed from OALIB's own .env files.
#
# The measured agreement the sweep asserts against is recorded in test/README.md -- update it there
# when a tolerance changes, so a loosened bound shows up in a diff.
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

            @testset "wavenumbers spanning several records" begin
                # LRecordLength is `MAX(2*Nfreq, 2*NzTab, 32, 3*nmedia)` longwords, and each complex
                # wavenumber takes two -- so few receivers plus many modes forces KRAKEN to split
                # the wavenumbers across records. Reading them all from one record, as both
                # published readers do, silently truncates here.
                result = KR.run_fortran_kraken(
                    UnderwaterEnv(munk_env()...), 200.0; rd=range(0.0, 5000.0; length=51), keep_files=true
                )
                try
                    lrecl = open(f -> Int(read(f, Int32)), joinpath(result.dir, "case.mod"))
                    records = 1 + (2 * result.nmodes - 1) ÷ lrecl
                    @test records > 1                       # the case is actually exercising it
                    @test result.nmodes > lrecl ÷ 2
                    @test length(result.kᵣ) == result.nmodes
                    # Wavenumbers descend across the whole set: a truncated read leaves zeros or
                    # garbage past the first record and breaks this.
                    @test issorted(real.(result.kᵣ); rev=true)
                    @test all(>(0), real.(result.kᵣ))
                    grp = result.grp
                    @test last(grp.m) > lrecl ÷ 2           # the .prt reaches past the first record
                    @test maximum(abs.(real.(result.kᵣ[grp.m]) .- real.(grp.kᵣ)) ./ abs.(real.(grp.kᵣ))) < 1e-6
                finally
                    rm(result.dir; recursive=true, force=true)
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

        @testset "comparison utility" begin
            pekeris = UnderwaterEnv(pekeris_env()...)

            @testset "canonical Pekeris agreement" begin
                c = KR.compare_with_fortran(pekeris, 100.0)
                @test c.n_julia == 5
                @test c.n_fortran == 5
                @test c.n_compared == 5
                @test KR.max_kr_reldiff(c) < 1e-5
                @test KR.min_mode_corr(c) > 0.999
                @test length(c.kr_absdiff) == length(c.kr_reldiff) == length(c.mode_corr) == 5
                @test c.kr_absdiff ≈ abs.(c.kr_julia .- c.kr_fortran)
                @test all(0 .<= c.mode_corr .<= 1 + 1e-12)
                @test c.group_speed_reldiff === nothing      # opt-in
                @test isempty(c.warnings)
            end

            @testset "summary table" begin
                text = sprint(show, MIME"text/plain"(), KR.compare_with_fortran(pekeris, 100.0))
                @test occursin("FortranComparison at 100.0 Hz", text)
                @test occursin("5 Julia modes vs 5 Fortran modes", text)
                @test occursin("max relative wavenumber difference", text)
                @test occursin("min mode-shape correlation", text)
                # One row per compared mode, plus header, title and the two summary lines.
                @test count(==('\n'), text) == 5 + 4
                # The compact form is what shows up inside containers.
                @test sprint(show, KR.compare_with_fortran(pekeris, 100.0)) ==
                    "FortranComparison(100.0 Hz, 5 vs 5 modes)"
            end

            @testset "nmodes caps the comparison" begin
                c = KR.compare_with_fortran(pekeris, 100.0; nmodes=3)
                @test c.n_julia == 5 && c.n_fortran == 5     # both sides still reported in full
                @test c.n_compared == 3
                @test length(c.kr_reldiff) == 3
            end

            @testset "mode correlation is sign- and scale-invariant" begin
                z = collect(0.0:1.0:100.0)
                a = sinpi.(z ./ 100)
                @test KR.mode_correlation(z, a, z, a) ≈ 1.0
                @test KR.mode_correlation(z, a, z, -7.5 .* a) ≈ 1.0      # sign and normalization
                @test KR.mode_correlation(z, a, z, sinpi.(2 .* z ./ 100)) < 1e-10   # orthogonal
                # Different grids: the Julia modes get resampled onto the Fortran depths.
                zf = collect(0.0:0.37:100.0)
                @test KR.mode_correlation(z, a, zf, sinpi.(zf ./ 100)) > 0.9999
            end

            @testset "group speeds against krakenc" begin
                # Off by default because it costs a ForwardDiff pass through the whole solver and
                # because the jll's kraken.exe reports no group speeds -- the comparison silently
                # re-runs with krakenc.exe to get them.
                c = KR.compare_with_fortran(pekeris, 100.0; group_speeds=true)
                @test c.group_speed_julia !== nothing
                @test c.group_speed_fortran !== nothing
                @test length(c.group_speed_reldiff) == 5
                @test all(1400.0 .< c.group_speed_julia .< 1600.0)
                # The .prt prints group speed with G14.6, and the highest-order mode sits nearest
                # the bottom cutoff where the two solvers' meshes disagree most; 1e-3 covers both.
                @test KR.max_group_speed_reldiff(c) < 1e-3
                @test occursin("max relative group-speed difference", sprint(show, MIME"text/plain"(), c))
            end

            @testset "a mode-count mismatch is reported, not thrown" begin
                # Kraken.jl searches kr in [ω/(0.9999 cb), max(ω/c)] while the .env asks KRAKEN for
                # phase speeds up to cb exactly, so Kraken.jl's window is very slightly narrower and
                # can miss a mode right at cutoff. That must surface as a number, not an exception.
                c = KR.compare_with_fortran(UnderwaterEnv(munk_env()...), 100.0)
                @test c.n_compared == min(c.n_julia, c.n_fortran)
                @test abs(c.n_julia - c.n_fortran) <= 1
                @test KR.max_kr_reldiff(c) < 1e-5
                @test KR.min_mode_corr(c) > 0.999
                if c.n_julia != c.n_fortran
                    @test occursin("mode-count mismatch", sprint(show, MIME"text/plain"(), c))
                end
            end
        end

        # --- cross-validation regression suite ---------------------------------------------
        #
        # Adding an environment is one row. `kr_rtol` and `corr_min` are per-environment because the
        # error is genuinely environment-dependent: a two-point isovelocity layer is nearly exact,
        # while a stack of gradient layers accumulates discretization error. Each bound is roughly
        # 3-10x the measured worst case (recorded in test/README.md) -- tight enough that a real
        # regression trips it, loose enough not to flake between binaries and platforms.
        #
        # `count_slack` is how many modes the two solvers may disagree on. It is 0 everywhere except
        # munk, where Kraken.jl finds one fewer mode at some frequencies: `bisection` searches phase
        # speeds up to `0.9999 cb` while the .env asks KRAKEN for up to `cb` exactly, so a mode
        # sitting right at the bottom cutoff can fall outside Kraken.jl's window.
        regression_cases = [
            (name="pekeris", build=() -> UnderwaterEnv(pekeris_env()...), kr_rtol=1e-6, corr_min=0.9999, count_slack=0),
            (
                name="one_layer",
                build=() -> UnderwaterEnv(one_layer_env()...),
                kr_rtol=1e-6,
                corr_min=0.9999,
                count_slack=0,
            ),
            (
                name="one_layer_slope",
                build=() -> UnderwaterEnv(one_layer_slope_env()...),
                kr_rtol=2e-5,
                corr_min=0.9999,
                count_slack=0,
            ),
            (
                name="two_layer_slope",
                build=() -> UnderwaterEnv(two_layer_slope_env()...),
                kr_rtol=1e-4,
                corr_min=0.9995,
                count_slack=0,
            ),
            (name="munk", build=() -> UnderwaterEnv(munk_env()...), kr_rtol=5e-5, corr_min=0.9999, count_slack=1),
        ]

        regression_freqs = [25.0, 50.0, 100.0, 200.0, 400.0]

        @testset "cross-validation against kraken.exe" begin
            @testset "$(case.name)" for case in regression_cases
                @testset "$(freq) Hz" for freq in regression_freqs
                    c = KR.compare_with_fortran(case.build(), freq)

                    # Both solvers must find something, or the comparison is vacuous.
                    @test c.n_julia > 0
                    @test c.n_fortran > 0
                    @test abs(c.n_julia - c.n_fortran) <= case.count_slack
                    @test c.n_compared == min(c.n_julia, c.n_fortran)

                    @test KR.max_kr_reldiff(c) < case.kr_rtol
                    @test KR.min_mode_corr(c) > case.corr_min

                    # Wavenumbers are ordered and physical: a trapped mode has phase speed between
                    # the slowest sound speed and the bottom half-space speed.
                    env = case.build()
                    @test issorted(c.kr_julia; rev=true)
                    @test all(2π * freq / env.cb .< c.kr_julia .< 2π * freq / minimum(env.c.c))

                    isempty(c.warnings) || @info "kraken.exe warnings for $(case.name) at $freq Hz" c.warnings
                end
            end
        end

        # --- AD gradients against Fortran (plan task 4.7) ------------------------------------
        #
        # Milestone 4's gradients are already checked against ForwardDiff and FiniteDiff in
        # test/reverse_ad_tests.jl. Those are strong tests of the *rules*, but every one of them
        # differentiates the same `det_sturm`: an error in the determinant itself moves the rrule,
        # ForwardDiff and FiniteDiff together and none of them notices. The two checks below use
        # kraken.exe as the oracle, so nothing they measure comes from Kraken.jl's differentiation.
        #
        # Fortran's mesh is pinned for both. KRAKEN's automatic mesh (NMESH = 0) is picked per
        # medium and is too coarse across two_layer_slope's 20 m layers to give an accurate
        # numerical dω/dkᵣ -- the group-speed disagreement there is 3.4e-3 on the automatic mesh
        # and 1.2e-4 on this one, a 29x improvement that tightening Kraken.jl's own tolerances
        # does not touch (measured 2026-08-09). Pinning also stops the mesh moving underneath a
        # finite difference, which would make the two perturbed runs incomparable.
        fine_mesh = 4000

        @testset "AD against Fortran" begin
            @testset "group speeds vs the .prt table" begin
                # KRAKEN computes group speeds numerically and prints them; Kraken.jl gets them by
                # differentiating the solver. Agreement is evidence about the derivative that does
                # not pass through any Kraken.jl AD code.
                #
                # Measured max relative difference on the pinned mesh (2026-08-09):
                #   pekeris 1.8e-6 | one_layer 2.6e-6 | one_layer_slope 3.0e-6
                #   two_layer_slope 1.2e-4 | munk 8.8e-5
                # The bound is the plan's 0.1%, which leaves 8-500x of headroom -- deliberately
                # loose, because this compares against whichever Fortran binary is present.
                @testset "$(case.name)" for case in regression_cases
                    c = KR.compare_with_fortran(case.build(), 100.0; group_speeds=true, nmesh=fine_mesh)

                    if c.group_speed_reldiff === nothing
                        # No binary here reports group speeds (AcousticsToolbox_jll's kraken.exe
                        # writes 0.00000 for every VG). compare_with_fortran already retries with
                        # krakenc; if that also fails there is nothing to compare against.
                        @info "no group speeds reported for $(case.name); skipping" maxlog = 5
                        @test_skip false
                    else
                        @test KR.max_group_speed_reldiff(c) < 1e-3
                        @test all(isfinite, c.group_speed_julia)
                        # A trapped mode's group speed is bounded by the medium sound speeds.
                        @test all(1000.0 .< c.group_speed_julia .< 2000.0)
                    end
                end
            end

            @testset "gradients vs finite-differenced kraken.exe" begin
                # The sharper of the two checks: perturb one environment parameter, write two .env
                # files, run the Fortran binary on each and central-difference its Re(kᵣ). That is
                # a gradient oracle for *any* parameter, not just frequency, and it is what caught
                # the layer-thickness derivatives being 15% off and sign-flipped during 4.3 while
                # ForwardDiff agreed with Zygote to 1e-9.

                # θ -> UnderwaterEnv, the same parameterizations test/reverse_ad_tests.jl uses.
                # Kept local rather than shared: these files are `include`d independently, and a
                # cross-file dependency would make either one unrunnable on its own.
                pekeris_θ(θ) = UnderwaterEnv(pekeris_env(; c0=θ[1], cb=θ[2], ρ0=θ[3], ρb=θ[4], depth=θ[5])...)
                onelayer_θ(θ) = UnderwaterEnv(
                    one_layer_env(; c0=θ[1], c1=θ[2], cb=θ[3], ρ0=θ[4], ρ1=θ[5], ρb=θ[6], h0=θ[7], h1=θ[8])...
                )

                # Tight tolerances so the comparison measures the rules and not the root solver --
                # see the TOL docstring in test/reverse_ad_tests.jl.
                tol = (abstol=1e-10, reltol=1e-10)
                kr_of(envf, θ, freq, mode) = kraken_jl(envf(θ), freq; tol...).kr[mode]

                # Re(kᵣ) straight out of kraken.exe, at the best precision the run reports: the
                # .prt's ten printed digits where it lists the mode, the .mod's single precision
                # otherwise.
                function fortran_kr(envf, θ, freq, mode)
                    ref = KR.run_fortran_kraken(envf(θ), freq; nmesh=fine_mesh)
                    return KR._best_fortran_kr(ref)[mode]
                end

                function fortran_derivative(envf, θ, k, freq, mode, h)
                    θp = copy(θ)
                    θp[k] *= (1 + h)
                    θm = copy(θ)
                    θm[k] *= (1 - h)
                    return (fortran_kr(envf, θp, freq, mode) - fortran_kr(envf, θm, freq, mode)) / (2 * θ[k] * h)
                end

                # Each row carries its own step, because the right step is set by the size of the
                # derivative rather than by the parameter. The .prt gives ten digits of kᵣ ~ 0.42,
                # so a difference below ~1e-10 is quantization noise: ∂kᵣ/∂h1 is 3.7e-8, and at
                # h = 1e-3 the two runs differ by only ~14 units in the last printed place, which
                # shows up as a 4.6% error. At h = 1e-2 the same row lands at 0.21%. Stepping the
                # other way is not free either -- cb at h = 1e-1 puts the half-space *below* the
                # water column and the run is rejected as unphysical before it starts.
                #
                # `rel` is the measured Zygote-vs-Fortran difference (2026-08-09); the assertion is
                # the plan's flat 1%, which every row clears by at least 4x.
                gradient_cases = [
                    (
                        env="pekeris",
                        envf=pekeris_θ,
                        θ=[1500.0, 1600.0, 1000.0, 1500.0, 100.0],
                        rows=[
                            (name="c0 (sound speed)", k=1, h=1e-3, rel=6.1e-7),
                            (name="ρ0 (density)", k=3, h=1e-2, rel=6.1e-5),
                            (name="depth (thickness)", k=5, h=1e-3, rel=1.8e-5),
                            (name="cb (control)", k=2, h=1e-3, rel=1.2e-4),
                        ],
                    ),
                    (
                        env="one_layer",
                        envf=onelayer_θ,
                        θ=[1500.0, 1550.0, 1600.0, 1000.0, 1500.0, 2000.0, 100.0, 20.0],
                        rows=[
                            (name="c1 (sound speed)", k=2, h=1e-3, rel=4.0e-4),
                            (name="ρ1 (density)", k=5, h=1e-3, rel=1.8e-4),
                            (name="h1 (thickness)", k=8, h=1e-2, rel=2.1e-3),
                            (name="c0 (control)", k=1, h=1e-3, rel=2.2e-6),
                        ],
                    ),
                ]

                @testset "$(case.env)" for case in gradient_cases
                    freq, mode = 100.0, 1
                    f = θ -> kr_of(case.envf, θ, freq, mode)

                    # One reverse-mode gradient gives every parameter at once; that is the whole
                    # point of the milestone, and it makes the sweep cheap.
                    g_reverse = Zygote.gradient(f, case.θ)[1]
                    g_forward = ForwardDiff.gradient(f, case.θ)

                    # Reverse mode must reproduce forward mode. Measured against the gradient's own
                    # scale, not entrywise: these gradients span four orders of magnitude (∂kᵣ/∂h1
                    # is 1.3e-4 of the largest entry), and an entrywise bound on the smallest
                    # components asks for agreement below the precision *either* method reaches --
                    # 1e-8 entrywise fails on h1 at 1.3e-7 while the two agree to 2.1e-11 on scale.
                    # Same reasoning as `relerr_norm` in test/reverse_ad_tests.jl.
                    # Measured 2026-08-09: pekeris 9.6e-12, one_layer 2.1e-11.
                    @test maximum(abs.(g_reverse .- g_forward)) / maximum(abs, g_forward) < 1e-9

                    @testset "$(row.name)" for row in case.rows
                        d_fortran = fortran_derivative(case.envf, case.θ, row.k, freq, mode, row.h)

                        # The independent assertion: this is the only check in the milestone that
                        # would survive Zygote and ForwardDiff being wrong in the same way, since
                        # `d_fortran` comes out of a separate binary. Loose by necessity -- see the
                        # per-row steps above for what the .prt's ten digits cost.
                        @test g_reverse[row.k] ≈ d_fortran rtol = 1e-2

                        # The derivative must be nonzero, or "agrees with Fortran" is vacuous.
                        @test abs(g_reverse[row.k]) > 0
                        # And it must have the sign Fortran gives it -- the failure mode 4.3 hit.
                        @test sign(g_reverse[row.k]) == sign(d_fortran)
                    end
                end
            end
        end

        # --- .env reader ---------------------------------------------------------------------

        @testset "env reader" begin
            # The checked-in files are the ones whose contents we know independently, because the
            # standard-environment builders produce the same environments. Parsing them and getting
            # the same cross-validation accuracy as the programmatic versions is the strongest
            # available evidence that the reader is faithful.
            reader_cases = [
                (file="Pekeris_AV", nmedia=1, depth=100.0, cb=1600.0, ρb=1500.0, nmodes=5),
                (file="onelayer_AV", nmedia=2, depth=120.0, cb=1600.0, ρb=2000.0, nmodes=5),
                (file="onelayer_slope_AV", nmedia=2, depth=120.0, cb=1600.0, ρb=2000.0, nmodes=5),
                (file="twolayer_slope_AV", nmedia=3, depth=140.0, cb=1800.0, ρb=2000.0, nmodes=10),
            ]

            @testset "$(case.file)" for case in reader_cases
                parsed = KR.read_env_file(joinpath(@__DIR__, "standard_envs", case.file * ".env"))
                @test length(parsed.env.layer_depth) == case.nmedia
                @test parsed.env.depth == case.depth
                @test parsed.env.cb == case.cb
                # Densities come back in kg/m³, not the file's g/cm³.
                @test parsed.env.ρb == case.ρb
                @test parsed.freqs == [100.0]
                @test parsed.clow == 1400.0

                c = KR.compare_with_fortran(parsed.env, parsed.freqs[1])
                @test c.n_julia == case.nmodes
                @test c.n_fortran == case.nmodes
                @test KR.max_kr_reldiff(c) < 1e-5
                @test KR.min_mode_corr(c) > 0.999
            end

            @testset "round trip through the writer" begin
                # write_env_file -> read_env_file must be the identity on the environment.
                for build in (pekeris_env, one_layer_env, one_layer_slope_env, two_layer_slope_env, munk_env)
                    original = UnderwaterEnv(build()...)
                    mktempdir() do dir
                        path = KR.write_env_file(joinpath(dir, "rt"), original, 100.0)
                        parsed = KR.read_env_file(path).env
                        @test parsed.layer_depth ≈ original.layer_depth
                        @test parsed.depth ≈ original.depth
                        @test parsed.cb ≈ original.cb
                        @test parsed.ρb ≈ original.ρb
                        @test parsed.c.c ≈ original.c.c
                        @test parsed.ρ.ρ ≈ original.ρ.ρ
                        @test -parsed.c.z ≈ -original.c.z
                    end
                end
            end

            @testset "unsupported features are named, not approximated" begin
                # ssp2.env declares n²-linear over a varying profile, which is genuinely a different
                # problem from the C-linear one Kraken.jl solves.
                err = try
                    KR.read_env_file(joinpath(@__DIR__, "standard_envs", "ssp2.env"))
                    nothing
                catch e
                    e
                end
                @test err isa KR.UnsupportedEnvFeature
                @test err.feature == "SSP interpolation"
                @test occursin("n²-linear", sprint(showerror, err))
                # ...but the environment is still recoverable for inspection.
                @test KR.read_env_file(joinpath(@__DIR__, "standard_envs", "ssp2.env"); strict=false).env isa
                    UnderwaterEnv

                # A cubic spline through a two-point medium *is* the straight line, so declaring 'S'
                # over the checked-in two-point files is not a reason to reject them.
                @test KR.read_env_file(joinpath(@__DIR__, "standard_envs", "Pekeris_AV.env")).topopt[1] == 'S'
            end

            @testset "a non-KRAKEN deck is rejected cleanly" begin
                mktempdir() do dir
                    path = joinpath(dir, "bellhop.env")
                    write(path, "'A ray deck'\n'NMNR'\n50.0\n-14.66, 14.66, 44\n")
                    err = try
                        KR.read_env_file(path)
                        nothing
                    catch e
                        e
                    end
                    @test err isa KR.MalformedEnvFile
                    @test occursin("NMEDIA", sprint(showerror, err))
                end
            end
        end

        # --- Acoustics Toolbox's own test cases ------------------------------------------------
        #
        # The toolbox tree is GPL-3 while this package is MIT, so its .env files are read in place
        # rather than vendored. Point KRAKEN_OALIB_TESTS at a checkout to run these; without one
        # they skip, which is what happens in CI.
        oalib_tree = get(ENV, "KRAKEN_OALIB_TESTS", "/Users/arielv/programs/AcousticsToolboxOALIB/tests")

        if isdir(oalib_tree)
            @testset "Acoustics Toolbox test cases" begin
                # A spread of shapes: a deep multi-layer Atlantic profile, a shallow penetrable
                # wedge slice, a 44-mode Pekeris waveguide, and the Munk profile at 102 modes.
                oalib_cases = [
                    (file="3DAtlantic/lante02.env", freq=50.0),
                    (file="3DAtlantic/lanta36.env", freq=50.0),
                    (file="Bellhop3DTests/PenetrableWedge/pwedge2d.env", freq=10.0),
                    (file="TLslices/pekeris.env", freq=10.0),
                    (file="Munk/MunkB_eigenray.env", freq=50.0),
                ]

                @testset "$(case.file)" for case in oalib_cases
                    path = joinpath(oalib_tree, case.file)
                    if !isfile(path)
                        @info "Not present in this Acoustics Toolbox checkout — skipped." path
                        continue
                    end
                    parsed = KR.read_env_file(path)
                    c = KR.compare_with_fortran(parsed.env, case.freq)
                    @test c.n_julia > 0
                    @test abs(c.n_julia - c.n_fortran) <= 1
                    @test KR.max_kr_reldiff(c) < 1e-4
                    @test KR.min_mode_corr(c) > 0.999
                end

                @testset "categorized report over the whole tree" begin
                    report = KR.categorize_env_tree(oalib_tree)
                    @test report.total > 100
                    @test length(report.supported) > 20
                    # Every rejection carries a reason; that list is the Milestone 5/6 backlog.
                    @test !isempty(report.unsupported)
                    @test all(!isempty, values(report.unsupported))
                    @test haskey(report.unsupported, "attenuation")
                    text = sprint(KR.print_env_tree_report, report)
                    @test occursin("Scanned $(report.total) .env files", text)
                    @test occursin("attenuation", text)
                    @info "Acoustics Toolbox coverage\n" * text
                end
            end
        else
            @info "No Acoustics Toolbox tree at $oalib_tree — its test cases are skipped. " *
                "Set KRAKEN_OALIB_TESTS to a checkout to run them."
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
