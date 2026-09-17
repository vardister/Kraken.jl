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

    # Reading a .prt does not need a Fortran binary, so this runs everywhere -- which matters,
    # because the format it pins is only ever produced by the Linux build.
    @testset "Group Speed table with an unfilled VG column" begin
        # Byte-for-byte the table AcousticsToolbox_jll's Linux kraken.exe writes for Pekeris_AV at
        # 100 Hz. Its group-speed array is uninitialised, and Fortran's F/G edit descriptors drop
        # the `E` from a three-digit exponent, so the garbage prints as `0.968142-315`. Parsing
        # that as end-of-table truncated the table to the single mode whose VG happened to be
        # clean, silently comparing mode 2 against mode 1's reference wavenumber.
        prt = """
 Nominal Frequency =   100.0    Hz

    I       k             alpha     Phase Speed     Group Speed
          (1/m)           (1/m)        (m/s)           (m/s)
    1  0.4179077218       0.0       1503.486291      0.968142-315
    2  0.4149628502       0.0       1514.156100       0.00000
    3  0.4099680381       0.0       1532.603697      0.106098-152
    4  0.4028501244       0.0       1559.683100      0.212024+162
    5  0.3938075024       0.0       1595.496599      0.234263-305
 _________________________________________________
"""
        mktempdir() do dir
            path = joinpath(dir, "garbage.prt")
            write(path, prt)
            grp = KR.read_grp(path)

            # The whole table, not just the rows that printed cleanly.
            @test grp.m == 1:5
            @test length(grp.kᵣ) == 5
            @test real(grp.kᵣ[1]) ≈ 0.4179077218 atol = 1e-9
            @test real(grp.kᵣ[5]) ≈ 0.3938075024 atol = 1e-9
            @test grp.phase_speed[1] ≈ 1503.486291 atol = 1e-6

            # Garbage is still garbage: the column must not be mistaken for real group speeds.
            @test !KR.has_group_speeds(grp)
        end
    end

    @testset "Fortran exponent spellings" begin
        @test KR._parse_fortran_float("1503.486291") ≈ 1503.486291
        @test KR._parse_fortran_float("-0.25") ≈ -0.25          # a leading sign is not an exponent
        @test KR._parse_fortran_float("0.968142-315") ≈ 0.968142e-315
        @test KR._parse_fortran_float("0.212024+162") ≈ 0.212024e162
        @test KR._parse_fortran_float("1.5D-3") ≈ 1.5e-3
        @test KR._parse_fortran_float("not a number") === nothing
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
            # `tempdir()` is shared with the rest of the machine, and not every entry in it is
            # readable: an Ubuntu CI runner has a root-only `/tmp/snap-private-tmp`, where `isfile`
            # raises EACCES rather than returning false. Anything we cannot stat is by definition
            # not one of our run directories, so treat a throw as "no".
            run_dirs() = filter(readdir(tempdir(); join=true)) do path
                try
                    isdir(path) && isfile(joinpath(path, "case.env"))
                catch
                    false
                end
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

        # --- attenuation against kraken.exe (plan task 5.3) ----------------------------------
        #
        # `Im(kᵣ)` needs tolerances of its own, and they are looser than the wavenumber ones by more
        # than the extra digits alone would explain. Three things stack up (see `max_alpha_reldiff`):
        # it is the small part of the complex number, the only Fortran source for it is the
        # single-precision `.mod`, and neither solver Richardson-extrapolates it — both evaluate the
        # perturbation on their own coarsest mesh.
        #
        # On top of that there is a real limit of the method, measured below and worth stating
        # plainly, because it sets what "agreement" can even mean here:
        #
        #   **Both solvers are first-order perturbation methods, and they agree to first order
        #   exactly** — the `α → 0` test below drives their disagreement down to 2.6e-6, which is
        #   the `.mod`'s single-precision floor. What they do *not* share is which second-order
        #   terms they keep, and the half-space term is where that shows: it goes as
        #   `Im √(γ² + 2iω α_b/c_b)`, whose expansion parameter is `2ω α_b/(c_b γ²)` — not small at
        #   all for a strongly attenuating bottom near cutoff. Measured on `SedAtten/calibK.env`
        #   (0.5 dB/λ half-space, 250 Hz) that parameter reaches 0.15 and the two solvers differ by
        #   ~1% on `Im(kᵣ)`; with the loss in the *water* instead, at the same 0.5 dB/λ, they agree
        #   to 9.3e-4. Neither is more right than the other: at that strength the first-order method
        #   itself is only good to about a percent, which is exactly why `krakenc.exe` — the full
        #   complex solve, plan task 5.6 — exists.
        #
        # So the plan's blanket 1e-3 target on the imaginary part is met for weakly attenuating
        # environments and is *not* achievable for strongly attenuating half-spaces by any
        # implementation of this method. The per-case bounds below are the usual 3-10x over the
        # measured worst case, and the measurements are recorded in test/README.md.
        # Built through `pekeris_env`'s own `α0`/`αb` keywords (task 5.5) rather than by poking the
        # matrices, so this is also a check that those keywords land in the columns they claim.
        function lossy_pekeris(; αb=0.0, αp=0.0, units=:dB_per_wavelength)
            return UnderwaterEnv(pekeris_env(; α0=αp, αb=αb)...; atten_units=units)
        end

        @testset "M5.3: modal attenuation against kraken.exe" begin
            # Measured on this waveguide at 100 Hz, 2026-08-09 (max relative difference over all
            # modes). Note how the attenuation column tracks `v = 2ω α_b/(c_b γ²)` and not the
            # attenuation itself: ten times the half-space loss costs a hundred times the agreement.
            #
            #   case                    Re(kᵣ)     Im(kᵣ)
            #   lossless control        1.7e-9     0, exactly, on both sides
            #   water 0.5 dB/λ          1.3e-4     4.0e-3
            #   half-space 0.05 dB/λ    1.3e-5     7.9e-4
            #   half-space 0.5 dB/λ     5.5e-4     1.0e-1
            atten_cases = [
                (name="lossless control", env=() -> lossy_pekeris(), kr_rtol=1e-6, α_rtol=1e-12, lossy=false),
                (name="water 0.5 dB/λ", env=() -> lossy_pekeris(; αp=0.5), kr_rtol=1e-3, α_rtol=2e-2, lossy=true),
                (
                    name="half-space 0.05 dB/λ",
                    env=() -> lossy_pekeris(; αb=0.05),
                    kr_rtol=1e-4,
                    α_rtol=5e-3,
                    lossy=true,
                ),
                (name="half-space 0.5 dB/λ", env=() -> lossy_pekeris(; αb=0.5), kr_rtol=3e-3, α_rtol=3e-1, lossy=true),
            ]

            @testset "$(case.name)" for case in atten_cases
                env = case.env()
                c = KR.compare_with_fortran(env, 100.0)

                @test c.n_julia == c.n_fortran
                @test c.n_julia > 0
                @test KR.max_kr_reldiff(c) < case.kr_rtol
                @test KR.min_mode_corr(c) > 0.999
                @test KR.max_alpha_reldiff(c) < case.α_rtol

                if case.lossy
                    # Both solvers must actually report loss, and with the same sign convention:
                    # a decaying mode has Im(kᵣ) < 0. Getting this backwards is the failure mode
                    # that a relative-difference test alone would not catch.
                    @test all(c.alpha_julia .< 0)
                    @test all(c.alpha_fortran .< 0)
                    @test eltype(kraken_jl(env, 100.0).kr) === ComplexF64
                else
                    # Fortran reports an exactly zero imaginary part for a lossless run, and so must
                    # Kraken.jl -- this is the control that says the machinery adds nothing.
                    @test all(iszero, c.alpha_fortran)
                    @test all(iszero, c.alpha_julia)
                    @test eltype(kraken_jl(env, 100.0).kr) === Float64
                end
            end
        end

        @testset "M5.5: the lossy standard environments cross-validate" begin
            # Task 5.5 put `α0`/`αb` keywords on `pekeris_env` and `α0`/`α1`/`αb` on `one_layer_env`.
            # Anything the package ships as a canned environment has to agree with kraken.exe, or it
            # is a trap for the first person who uses it — so the lossy variants join the regression
            # list rather than living only in the docs.
            #
            # `one_layer_env(; α1=…)` is the important row: it is the only case here with a lossy
            # *interior* medium, and it is what exposed the discontinuous-quadrature error that task
            # 5.7 then fixed. Measured before and after that fix, so the effect is on the record:
            #
            #                                          Re(kᵣ)    Im(kᵣ) before   Im(kᵣ) after
            #   pekeris_env(α0=0.2)                    2.1e-5    7.9e-4          8.1e-4
            #   pekeris_env(αb=0.2)                    1.6e-4    5.8e-2          5.8e-2
            #   one_layer_env(α1=0.4)                  3.0e-5    8.1e-2          2.9e-3
            #   one_layer_env(α0=.1, α1=.4, αb=.1)     4.3e-5    8.7e-3          7.5e-4
            #
            # Exactly the pattern the fix predicts. The two `pekeris` rows are a single water medium
            # over a half-space, so they have no interior interface to straddle and are untouched —
            # their residual is the bottom-cutoff limit, which is a property of first-order
            # perturbation theory and not of the quadrature. The two `one_layer` rows have a lossy
            # sediment between two interfaces and improve 27x and 12x.
            std_lossy = [
                (
                    name="pekeris_env(α0=0.2)",
                    env=() -> UnderwaterEnv(pekeris_env(; α0=0.2)...),
                    kr_rtol=1e-4,
                    α_rtol=5e-3,
                ),
                (
                    name="pekeris_env(αb=0.2)",
                    env=() -> UnderwaterEnv(pekeris_env(; αb=0.2)...),
                    kr_rtol=1e-3,
                    α_rtol=2e-1,
                ),
                (
                    name="one_layer_env(α1=0.4)",
                    env=() -> UnderwaterEnv(one_layer_env(; α1=0.4)...),
                    kr_rtol=1e-4,
                    α_rtol=1e-2,
                ),
                (
                    name="one_layer_env(α0=0.1, α1=0.4, αb=0.1)",
                    env=() -> UnderwaterEnv(one_layer_env(; α0=0.1, α1=0.4, αb=0.1)...),
                    kr_rtol=1e-4,
                    α_rtol=5e-3,
                ),
            ]

            @testset "$(case.name)" for case in std_lossy
                env = case.env()
                @test is_lossy(env)
                c = KR.compare_with_fortran(env, 100.0)
                @test c.n_julia == c.n_fortran
                @test c.n_julia > 0
                @test KR.max_kr_reldiff(c) < case.kr_rtol
                @test KR.min_mode_corr(c) > 0.999
                @test KR.max_alpha_reldiff(c) < case.α_rtol
                @test all(c.alpha_julia .< 0)
                @test all(c.alpha_fortran .< 0)
            end

            # The keywords default to zero, so the canned environments are lossless unless asked --
            # this is what keeps every pre-Milestone-5 result in this file unchanged.
            for build in (pekeris_env, one_layer_env, one_layer_slope_env, two_layer_slope_env, munk_env)
                @test !is_lossy(UnderwaterEnv(build()...))
            end
        end

        @testset "M5.7: the perturbation integral converges at second order" begin
            # The tolerance above says the lossy-layer case is now accurate; this says *why*, which
            # is the part that will still be true after someone changes the mesh defaults.
            #
            # Before 5.7 a single trapezoid ran across the α jump at the top of the sediment, and the
            # error fell as O(h): 8.1e-2, 4.0e-2, 1.9e-2, 8.8e-3, 3.6e-3 at 20/40/80/160/320 points
            # per wavelength — a measured order of 1.12. Integrating medium by medium removes the
            # straddled interval, and the error should now fall as O(h²).
            env = UnderwaterEnv(one_layer_env(; α1=0.4)...)
            freq = 100.0

            function alpha_at(npw)
                props = AcousticProblemProperties(env, freq; n_per_wavelength=npw)
                cache = AcousticProblemCache(env, props)
                krc, ψ = inverse_iteration(find_kr(env, props, cache), env, props, cache; reltol=1e-10)
                return map(eachindex(krc)) do m
                    δ = Kraken.modal_attenuation(view(ψ, :, m), krc[m], env, props)
                    return imag(sqrt(complex(krc[m]^2) + δ))
                end
            end

            # The reference is Kraken.jl's *own* finely resolved answer, not kraken.exe's. That is
            # deliberate and it is the only way this measures what it claims to: `kraken.exe` runs on
            # its own automatic mesh and carries about 1.6e-3 of discretization error on this case,
            # which is a floor, not a slope. Measured against it the observed order flattens to 0.57,
            # 0.19, 0.05 as our error drops below Fortran's — an artifact of the reference, and
            # exactly the trap this comment exists to stop someone falling into.
            reference = alpha_at(320)[1]
            errs = [abs(alpha_at(npw)[1] - reference) / abs(reference) for npw in (20, 40, 80)]

            @test all(errs .> 0)
            @test issorted(errs; rev=true)

            # Order = log2 of the ratio between successive mesh doublings. Before 5.7 this was 1.12
            # (a single trapezoid straddling the α jump at the sediment top); after, it is 2.00.
            orders = [log2(errs[i] / errs[i + 1]) for i in 1:(length(errs) - 1)]
            @test all(orders .> 1.8)

            # ...and the agreement with kraken.exe at the *default* mesh, which is the number a user
            # actually meets. 8.1e-2 before 5.7, 2.9e-3 after.
            fortran = imag.(KR.run_fortran_kraken(env, freq).kᵣ)
            @test abs(alpha_at(20)[1] - fortran[1]) / abs(fortran[1]) < 1e-2
        end

        @testset "M5.3: the two solvers agree to first order in α" begin
            # The sharpest statement available about the attenuation, and the one that shows the
            # ~1% seen above is a second-order artifact rather than a defect: shrink the
            # attenuation and the disagreement shrinks *faster* than linearly, bottoming out at the
            # single precision of the .mod file. A formula error would leave a floor proportional
            # to α instead.
            errs = map((0.5, 0.05, 0.005)) do αb
                c = KR.compare_with_fortran(lossy_pekeris(; αb=αb), 250.0)
                return c.alpha_reldiff[1]
            end

            @test errs[3] < errs[2] < errs[1]        # monotone in α
            @test errs[1] / errs[2] > 10             # faster than linear: it is the α² term
            @test errs[3] < 1e-4                     # and it lands at the .mod's precision floor
        end

        @testset "M5.3: the writer round-trips attenuation and its units" begin
            # A `.env` records attenuation *values* in column 6 and their *units* in column 3 of the
            # top-option string. Writing one without the other silently rescales every attenuation
            # in the file, so the two are tested together.
            @testset "$units" for units in
                                  (:nepers_per_m, :dB_per_m, :dB_per_kmHz, :dB_per_wavelength, :Q, :loss_parameter)
                env = lossy_pekeris(; αb=0.5, αp=0.02, units=units)
                text = KR.env_file_string(env, 100.0)
                char = only(filter(p -> p.second === units, collect(ATTENUATION_UNIT_CHARS))).first
                @test occursin("'CV$char'", text)

                mktempdir() do dir
                    path = KR.write_env_file(joinpath(dir, "lossy"), env, 100.0)
                    back = KR.read_env_file(path)
                    @test back.atten_units === units
                    @test back.env.αb ≈ env.αb
                    @test back.env.α.α ≈ env.α.α
                    @test back.env.cb ≈ env.cb
                end
            end

            # An `UnderwaterEnvFORTRAN` carries no units of its own, so the keyword supplies them.
            ssp, layers, sspHS = pekeris_env()
            sspHS[2, 5] = 0.3
            envf = UnderwaterEnvFORTRAN(ssp, layers, sspHS)
            @test occursin("'CVN'", KR.env_file_string(envf, 100.0; atten_units=:nepers_per_m))
            @test occursin("'CVW'", KR.env_file_string(envf, 100.0))

            # An explicitly supplied non-default units letter is the caller's choice and is kept.
            @test occursin("'CVQ'", KR.env_file_string(lossy_pekeris(; αb=0.5), 100.0; topopt="CVQ"))
        end

        @testset "M5.3: kraken.exe accepts the lossy files the writer emits" begin
            # Same check as the lossless writer suite at the bottom of this file: kraken.exe exits 0
            # even on a fatal error, so the .prt is the only honest signal. A wrong units letter or
            # a malformed attenuation column shows up here as an ERROR line.
            @testset "$units" for units in (:nepers_per_m, :dB_per_wavelength, :dB_per_kmHz)
                mktempdir() do dir
                    env = lossy_pekeris(; αb=units === :nepers_per_m ? 0.001 : 0.3, units=units)
                    KR.write_env_file(joinpath(dir, "lossy"), env, 100.0)
                    cmd = Cmd(`$(KR.kraken_cmd()) lossy`; dir=dir)
                    run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull))
                    report = read(joinpath(dir, "lossy.prt"), String)
                    @test isempty(filter(l -> occursin("ERROR", uppercase(l)), split(report, '\n')))
                    # The .prt echoes the units it decided on, which is the strongest confirmation
                    # that the letter we wrote is the one it read.
                    @test occursin("Attenuation units", report)
                end
            end
        end

        # Every boundary-condition pair and interpolation Kraken.jl solves (plan task 6.4), on the
        # Munk profile: one medium, so a cubic spline is allowed, and curved, so the three
        # interpolators genuinely differ. `topopt`/`botopt` are the letters `ReadTopOpt` and `TopBot`
        # have to see.
        option_cases = [
            (top=PressureRelease(), bottom=AcousticHalfspace(), interp=:c_linear, topopt="CVW", botopt=('A')),
            (top=RigidBoundary(), bottom=AcousticHalfspace(), interp=:c_linear, topopt="CRW", botopt=('A')),
            (top=PressureRelease(), bottom=RigidBoundary(), interp=:c_linear, topopt="CVW", botopt=('R')),
            (top=PressureRelease(), bottom=PressureRelease(), interp=:c_linear, topopt="CVW", botopt=('V')),
            (top=RigidBoundary(), bottom=RigidBoundary(), interp=:c_linear, topopt="CRW", botopt=('R')),
            (top=RigidBoundary(), bottom=PressureRelease(), interp=:n2_linear, topopt="NRW", botopt=('V')),
            (top=PressureRelease(), bottom=AcousticHalfspace(), interp=:n2_linear, topopt="NVW", botopt=('A')),
            (top=PressureRelease(), bottom=AcousticHalfspace(), interp=:cubic_spline, topopt="SVW", botopt=('A')),
        ]
        option_env(case) = UnderwaterEnv(munk_env()...; ssp_interp=case.interp, top_bc=case.top, bottom_bc=case.bottom)

        @testset "M6.4: the writer declares interpolation and boundary conditions" begin
            @testset "$(case.topopt) over $(case.botopt)" for case in option_cases
                env = option_env(case)
                text = KR.env_file_string(env, 10.0)
                @test occursin("'$(case.topopt)'", text)
                @test occursin("'$(case.botopt)' ", text)
                # Only a half-space has a record of its own, and only a half-space bounds the
                # trapped spectrum; a perfect bottom traps modes all the way down to kᵣ = 0.
                @test occursin("! ZB  CPB", text) == (case.botopt == 'A')
                @test occursin(KR._fmt(KR.PERFECT_BOTTOM_CHIGH), text) == (case.botopt != 'A')

                mktempdir() do dir
                    back = KR.read_env_file(KR.write_env_file(joinpath(dir, "options"), env, 10.0))
                    @test back.env.c.mode === case.interp
                    @test back.env.top_bc == case.top
                    @test back.env.bottom_bc == case.bottom
                    @test back.env.c.c ≈ env.c.c
                    @test -back.env.c.z ≈ -env.c.z
                    @test back.nmesh == [0]
                    case.botopt == 'A' && @test back.env.cb ≈ env.cb
                end
            end

            # An explicit letter is the caller's choice and survives, in every column.
            env = option_env(option_cases[end])
            @test occursin("'NVW'", KR.env_file_string(env, 10.0; topopt="NVW"))
            @test occursin("'SRQ'", KR.env_file_string(option_env(option_cases[1]), 10.0; topopt="SRQ"))
            text = KR.env_file_string(option_env(option_cases[3]), 10.0; botopt="A")
            @test occursin("'A' ", text) && occursin("! ZB  CPB", text)

            # An `UnderwaterEnvFORTRAN` says nothing about either, so the defaults stand.
            envf = UnderwaterEnvFORTRAN(munk_env()...)
            text = KR.env_file_string(envf, 10.0)
            @test occursin("'CVW'", text)
            @test occursin("'A' ", text)
        end

        @testset "M6.4: kraken.exe accepts every option the writer emits" begin
            # The .prt echoes each choice back, which confirms the letter we wrote is the one it
            # read -- not merely that the file parsed.
            echo = Dict(
                :c_linear => "C-Linear approximation",
                :n2_linear => "N2-Linear approximation",
                :cubic_spline => "Spline approximation",
                'A' => "ACOUSTO-ELASTIC half-space",
                'V' => "VACUUM",
                'R' => "Perfectly RIGID",
            )
            @testset "$(case.topopt) over $(case.botopt)" for case in option_cases
                ref = KR.run_fortran_kraken(option_env(case), 10.0; keep_files=true)
                try
                    @test ref.nmodes > 0
                    report = read(joinpath(ref.dir, "case.prt"), String)
                    @test occursin(echo[case.interp], report)
                    @test occursin(echo[case.botopt], report)
                    @test occursin(echo[case.topopt[2]], report)
                finally
                    rm(ref.dir; recursive=true, force=true)
                end
            end
        end

        # --- Milestone 6 options against kraken.exe (plan task 6.5) ---------------------------
        #
        # Measured 2026-09-16; the tables are in test/README.md under "Boundary conditions and SSP
        # interpolation validated against Fortran".

        @testset "M6.5: every boundary condition against kraken.exe" begin
            # Pekeris' isovelocity column isolates the boundary rows: nothing else about the problem
            # differs between the cases. 100 Hz is away from every mode cutoff (see the broken tests
            # below for why that matters).
            bc_cases = [
                (top=RigidBoundary(), bottom=AcousticHalfspace(), freq=100.0, nmodes=5, kr_rtol=1e-8),
                (top=PressureRelease(), bottom=RigidBoundary(), freq=100.0, nmodes=13, kr_rtol=1e-8),
                (top=PressureRelease(), bottom=PressureRelease(), freq=100.0, nmodes=13, kr_rtol=1e-7),
                (top=RigidBoundary(), bottom=RigidBoundary(), freq=100.0, nmodes=14, kr_rtol=1e-8),
                (top=RigidBoundary(), bottom=PressureRelease(), freq=100.0, nmodes=13, kr_rtol=1e-8),
                (top=PressureRelease(), bottom=PressureRelease(), freq=50.0, nmodes=6, kr_rtol=1e-8),
                # Mode 7 is grazing, kᵣ = 0.046 against ω/c = 0.209, so a small absolute error in kᵣ²
                # reads large relative to kᵣ: 6.1e-5 here, from 1.4e-10 at 100 Hz.
                (top=PressureRelease(), bottom=RigidBoundary(), freq=50.0, nmodes=7, kr_rtol=2e-4),
            ]
            @testset "$(case.top) over $(case.bottom), $(case.freq) Hz" for case in bc_cases
                env = UnderwaterEnv(pekeris_env()...; top_bc=case.top, bottom_bc=case.bottom)
                c = KR.compare_with_fortran(env, case.freq)
                @test c.n_julia == case.nmodes
                @test c.n_fortran == case.nmodes
                @test KR.max_kr_reldiff(c) < case.kr_rtol
                @test KR.min_mode_corr(c) > 0.9999
            end

            # A curved profile over a perfect bottom, so the boundary rows meet a varying column.
            c = KR.compare_with_fortran(
                UnderwaterEnv(munk_env()...; ssp_interp=:n2_linear, bottom_bc=RigidBoundary()), 10.0
            )
            @test c.n_julia == c.n_fortran == 66
            @test KR.max_kr_reldiff(c) < 5e-5
            @test KR.min_mode_corr(c) > 0.9999

            # A perfect bottom can crash the solve when a mode sits at kᵣ = 0 -- for this 100 m
            # column the cutoffs are the multiples of c/2D = 7.5 Hz (vacuum) and the odd multiples
            # of 3.75 Hz (rigid). `richard_extrap` takes the square root of an extrapolated kᵣ² that
            # has crossed zero; KRAKEN discards such a mode instead. Measured over 20-200 Hz in
            # 0.25 Hz steps: 10 of 721 frequencies fail for vacuum, 6 for rigid, three of those a
            # SingularException in inverse iteration at 183.25-183.75 Hz. Every failure is at or
            # beside a cutoff, but not every cutoff fails. Plan task 6.7.
            solves(env, freq) =
                try
                    kraken_jl(env, freq)
                    true
                catch
                    false
                end
            @test_broken solves(UnderwaterEnv(pekeris_env()...; bottom_bc=PressureRelease()), 75.0)
            @test_broken solves(UnderwaterEnv(pekeris_env()...; bottom_bc=RigidBoundary()), 93.75)
            @test_broken solves(UnderwaterEnv(pekeris_env()...; bottom_bc=RigidBoundary()), 183.5)
        end

        @testset "M6.5: every SSP interpolation against kraken.exe" begin
            # A coarse, strongly curved duct: five samples 50 m apart, so the three interpolants
            # genuinely disagree between them. On `munk_env`'s 100 m sampling of a smooth profile
            # they do not -- at 25 Hz the n²-linear solve matched Fortran's n²-linear run to 6.4e-6
            # and Fortran's *C-linear* run to 6.8e-6, which cannot tell the options apart.
            duct_z = [0.0, 50.0, 100.0, 150.0, 200.0]
            duct_c = [1540.0, 1500.0, 1480.0, 1500.0, 1530.0]
            function duct(z=duct_z, c=duct_c; kw...)
                n = length(z)
                ssp = hcat(z, c, zeros(n), fill(1000.0, n), zeros(n), zeros(n))
                sspHS = [0.0 343.0 0.0 0.00121 0.0 0.0; 200.0 1700.0 0.0 1800.0 0.0 0.0]
                return UnderwaterEnv(ssp, [0.0 0.0 200.0], sspHS; kw...)
            end
            function reldiff(a, b)
                n = min(length(a), length(b))
                return maximum(abs.(a[1:n] .- b[1:n]) ./ abs.(b[1:n]))
            end

            interps = (:c_linear, :n2_linear, :cubic_spline)
            @testset "$(freq) Hz" for freq in (50.0, 100.0)
                fortran = Dict(
                    m => KR._best_fortran_kr(KR.run_fortran_kraken(duct(; ssp_interp=m), freq)) for m in interps
                )
                julia = Dict(m => Float64.(real.(kraken_jl(duct(; ssp_interp=m), freq).kr)) for m in interps)
                for m in interps
                    @test length(julia[m]) == length(fortran[m])
                end
                # The matrix of max rel Δkᵣ, Julia interpolation × Fortran interpolation. Measured:
                #                 C        N        S       (50 Hz; 100 Hz is within 4x of each entry)
                #   C-linear   7.7e-8   1.2e-4   1.8e-3
                #   n²-linear  1.2e-4   7.0e-8   1.8e-3
                #   spline     1.8e-3   1.7e-3   3.1e-4
                # C and N are told apart by three orders of magnitude, so a wrong interpolator fails.
                for m in (:c_linear, :n2_linear)
                    @test reldiff(julia[m], fortran[m]) < 1e-6
                    @test all(reldiff(julia[m], fortran[o]) > 1e-5 for o in interps if o !== m)
                end
                # The spline agrees only to 3e-4, and the reason is the end condition, not the solver
                # (checked below): Kraken.jl's spline is DataInterpolations' natural spline, KRAKEN's
                # `CSPLINE` is called with IBCBEG = IBCEND = 0, which `splinec.f90` documents as
                # not-a-knot. Plan task 6.8.
                @test reldiff(julia[:cubic_spline], fortran[:cubic_spline]) < 5e-4
                @test_broken reldiff(julia[:cubic_spline], fortran[:cubic_spline]) < 1e-6
                @test all(reldiff(julia[:cubic_spline], fortran[o]) > 1e-3 for o in (:c_linear, :n2_linear))
            end

            @testset "the spline gap is the end condition" begin
                # Sample each end condition's spline through the duct at 0.5 m and solve that as a
                # C-linear profile. The not-a-knot one reproduces Fortran's spline run; the natural
                # one reproduces Kraken.jl's. Measured at 50 Hz: 2.1e-7 and 2.0e-7, against 3.1e-4
                # crosswise.
                function spline_moments(x, y, not_a_knot)
                    n = length(x)
                    h = diff(x)
                    A = zeros(n, n)
                    b = zeros(n)
                    for i in 2:(n - 1)
                        A[i, (i - 1):(i + 1)] = [h[i - 1], 2 * (h[i - 1] + h[i]), h[i]]
                        b[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
                    end
                    if not_a_knot   # third derivative continuous across x[2] and x[n-1]
                        A[1, 1:3] = [-1 / h[1], 1 / h[1] + 1 / h[2], -1 / h[2]]
                        A[n, (n - 2):n] = [-1 / h[n - 2], 1 / h[n - 2] + 1 / h[n - 1], -1 / h[n - 1]]
                    else            # natural: zero second derivative at both ends
                        A[1, 1] = A[n, n] = 1
                    end
                    return A \ b
                end
                function spline_at(x, y, M, xq)
                    i = clamp(searchsortedlast(x, xq), 1, length(x) - 1)
                    h, a, b = x[i + 1] - x[i], x[i + 1] - xq, xq - x[i]
                    return M[i] * a^3 / (6h) +
                           M[i + 1] * b^3 / (6h) +
                           (y[i] / h - M[i] * h / 6) * a +
                           (y[i + 1] / h - M[i + 1] * h / 6) * b
                end
                dense = collect(0.0:0.5:200.0)
                sampled(not_a_knot) =
                    let M = spline_moments(duct_z, duct_c, not_a_knot)
                        Float64.(
                            real.(kraken_jl(duct(dense, [spline_at(duct_z, duct_c, M, z) for z in dense]), 50.0).kr)
                        )
                    end
                fortran_spline = KR._best_fortran_kr(KR.run_fortran_kraken(duct(; ssp_interp=:cubic_spline), 50.0))
                julia_spline = Float64.(real.(kraken_jl(duct(; ssp_interp=:cubic_spline), 50.0).kr))
                @test reldiff(sampled(true), fortran_spline) < 1e-6
                @test reldiff(sampled(false), julia_spline) < 1e-6
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

            @testset "M6.4: interpolation and boundary options are read" begin
                # ssp2.env declares n²-linear over a varying three-medium profile. Until Milestone 6
                # that was its refusal; now the interpolation is read, and what still stops it is
                # physics: its half-space is exactly as fast as its fastest water (1600.33 m/s), so
                # nothing is trapped.
                ssp2_path = joinpath(@__DIR__, "standard_envs", "ssp2.env")
                err = try
                    KR.read_env_file(ssp2_path)
                    nothing
                catch e
                    e
                end
                @test err isa KR.UnsupportedEnvFeature
                @test err.feature == "bottom half-space is not the fastest medium"
                ssp2 = KR.read_env_file(ssp2_path; strict=false)
                @test ssp2.env.c.mode === :n2_linear
                @test ssp2.nmesh == [0, 0, 0]

                # A cubic spline through two-point media *is* the straight line, so the checked-in
                # files that declare 'S' read as the C-linear problem they describe.
                pekeris = KR.read_env_file(joinpath(@__DIR__, "standard_envs", "Pekeris_AV.env"))
                @test pekeris.topopt[1] == 'S'
                @test pekeris.env.c.mode === :c_linear

                # A minimal deck with a choice of interpolation, top and bottom. `media` is a list of
                # (bottom depth, [(z, c), ...]) and the half-space record is written only for 'A'.
                function deck(topopt, botopt; media=[(100.0, [(0.0, 1500.0), (100.0, 1500.0)])], nmesh=0)
                    io = IOBuffer()
                    println(io, "'options probe'\n50.0\n$(length(media))\n'$topopt'")
                    for (bottom, rows) in media
                        println(io, "$nmesh  0.0  $bottom")
                        for (z, c) in rows
                            println(io, "  $z  $c  0.0  1.0  0.0  0.0 /")
                        end
                    end
                    println(io, "'$botopt'  0.0")
                    botopt == "A" && println(io, "  $(first(last(media)))  1600.0  0.0  1.5  0.0  0.0 /")
                    println(io, "1400.0  $(botopt == "A" ? 1600.0 : 1.0e7)\n10.0")
                    return String(take!(io))
                end
                read_deck(text) = mktempdir() do dir
                    path = joinpath(dir, "probe.env")
                    write(path, text)
                    KR.read_env_file(path)
                end
                refusal(text) =
                    try
                        read_deck(text)
                        nothing
                    catch e
                        e
                    end

                curved = [(100.0, [(0.0, 1500.0), (50.0, 1490.0), (100.0, 1510.0)])]
                two_curved = [curved[1], (200.0, [(100.0, 1520.0), (150.0, 1530.0), (200.0, 1560.0)])]
                two_linear = [(100.0, [(0.0, 1500.0), (100.0, 1490.0)]), (200.0, [(100.0, 1520.0), (200.0, 1560.0)])]

                for (char, mode) in (('C', :c_linear), ('N', :n2_linear), ('S', :cubic_spline))
                    @test read_deck(deck("$(char)VW", "A"; media=curved)).env.c.mode === mode
                end
                @test read_deck(deck("SVW", "A"; media=two_linear)).env.c.mode === :c_linear
                @test read_deck(deck("PVW", "A"; media=two_linear)).env.c.mode === :c_linear
                @test read_deck(deck("NVW", "A"; media=two_curved)).env.c.mode === :n2_linear

                for (char, bc) in (('V', PressureRelease()), ('R', RigidBoundary()))
                    @test read_deck(deck("C$(char)W", "A")).env.top_bc == bc
                end
                for (char, bc) in (('A', AcousticHalfspace()), ('V', PressureRelease()), ('R', RigidBoundary()))
                    parsed = read_deck(deck("CVW", string(char)))
                    @test parsed.env.bottom_bc == bc
                    @test parsed.chigh == (char == 'A' ? 1600.0 : 1.0e7)
                    # A perfect bottom has no record, so there is nothing to fill `sspHS` with.
                    char == 'A' || @test parsed.sspHS[2, 2:end] == zeros(5)
                end
                # A perfect bottom traps every mode, so it is not held to "fastest medium".
                @test read_deck(deck("CRW", "R"; media=curved)).env.bottom_bc == RigidBoundary()
                @test read_deck(deck("CVW", "A"; media=curved, nmesh=250)).nmesh == [250]

                # What is still refused is named, with the reason.
                err = refusal(deck("SVW", "A"; media=two_curved))
                @test err isa KR.UnsupportedEnvFeature
                @test err.feature == "SSP interpolation"
                @test occursin("cubic spline over 2 media", sprint(showerror, err))
                # ...but the environment is still recoverable for inspection, as the linear profile.
                mktempdir() do dir
                    path = joinpath(dir, "spline.env")
                    write(path, deck("SVW", "A"; media=two_curved))
                    @test KR.read_env_file(path; strict=false).env.c.mode === :c_linear
                end

                err = refusal(deck("PVW", "A"; media=curved))
                @test err isa KR.UnsupportedEnvFeature
                @test occursin("PCHIP", sprint(showerror, err))

                err = refusal(deck("CFW", "A"))
                @test err isa KR.UnsupportedEnvFeature
                @test err.feature == "top boundary"

                err = refusal(deck("CVW", "F"))
                @test err isa KR.UnsupportedEnvFeature
                @test err.feature == "bottom boundary"
                @test occursin("reflection-coefficient file", sprint(showerror, err))
            end

            @testset "M5.1: attenuation units are read off the top-option string" begin
                # `read_env_file` builds a minimal deck for each unit character so the mapping is
                # tested without needing a toolbox checkout. Column 3 of TOPOPT is the unit; the
                # value in SSP column 5 is carried through unchanged, because converting it needs a
                # frequency and the environment does not have one.
                function deck(topopt, αp, αb)
                    return """
                    'units probe'
                    100.0
                    1
                    '$topopt'
                    0  0.0  100.0
                         0.0  1500.0  0.0  1.0  $αp  0.0 /
                       100.0  1500.0  0.0  1.0  $αp  0.0 /
                    'A'  0.0
                       100.0  1600.0  0.0  1.5  $αb  0.0 /
                    1400.0  1600.0
                    10.0
                    """
                end
                mktempdir() do dir
                    for (char, units) in ATTENUATION_UNIT_CHARS
                        path = joinpath(dir, "units_$char.env")
                        write(path, deck("CV$char", 0.02, 0.5))
                        parsed = KR.read_env_file(path)
                        @test parsed.atten_units === units
                        @test parsed.env.atten_units === units
                        @test parsed.env.αb == 0.5
                        @test all(≈(0.02), parsed.env.α.α)
                        @test is_lossy(parsed.env)
                    end

                    # A zero attenuation is still lossless whatever the declared units are.
                    path = joinpath(dir, "lossless.env")
                    write(path, deck("CVW", 0.0, 0.0))
                    @test !is_lossy(KR.read_env_file(path).env)

                    # An unusable units character is rejected rather than defaulted -- `ReadTopOpt`
                    # calls ERROUT on exactly the same input, so a file we accepted here would be one
                    # kraken.exe refuses to run.
                    for bad in ("CVX", "CV")
                        path = joinpath(dir, "bad.env")
                        write(path, deck(bad, 0.0, 0.0))
                        err = try
                            KR.read_env_file(path)
                            nothing
                        catch e
                            e
                        end
                        @test err isa KR.MalformedEnvFile
                        @test occursin("attenuation-units", sprint(showerror, err))
                    end

                    # 'm' is the seventh convention -- dB/m with a power law -- and needs per-medium
                    # parameters this reader does not model. It must be named, not read as 'M'.
                    path = joinpath(dir, "powerlaw.env")
                    write(path, deck("CVm", 0.1, 0.1))
                    err = try
                        KR.read_env_file(path)
                        nothing
                    catch e
                        e
                    end
                    @test err isa KR.UnsupportedEnvFeature
                    @test err.feature == "power-law attenuation"
                end
            end

            @testset "M5.1: the bottom-option record ends after SIGMA" begin
                # `READ( ENVFile, * ) BotOpt, Sigma` takes two items and discards the rest of the
                # record. SedAtten/calibS_0.6dB.env writes `'A'  0.0 2.5 2000`; reading the *last*
                # number as the roughness made that file look like it had 2 km of interfacial
                # roughness and rejected it.
                mktempdir() do dir
                    path = joinpath(dir, "trailing.env")
                    write(
                        path,
                        """
                        'trailing values on the bottom record'
                        250.0
                        1
                        'CVW'
                        0  0.0  100.0
                             0.0  1500.0  0.0  1.0  0.0  0.0 /
                           100.0  1500.0  0.0  1.0  0.0  0.0 /
                        'A'  0.0 2.5 2000
                           100.0  1590.0  0.0  1.2  0.5  0.0 /
                        1400.0  1590.0
                        30.0
                        """,
                    )
                    parsed = KR.read_env_file(path)
                    @test parsed.sigmas == [0.0, 0.0]
                    @test parsed.env.αb == 0.5
                    @test parsed.env.cb == 1590.0
                end
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

                @testset "M5.3: SedAtten and the attenuation TL slice against kraken.exe" begin
                    # The milestone's named cases. `VolAtt` is deliberately absent and the reason is
                    # worth recording: every file in it declares an acousto-elastic half-space
                    # *above* the surface (`TopOpt(2:2) == 'A'`) and a bottom whose sound speed
                    # equals the water's, so there is no trapped spectrum to compare -- they are
                    # free-space transmission-loss cases, not modal ones. Two of them additionally
                    # ask for Thorp or Francois-Garrison volume attenuation, which is `TopOpt(4:4)`
                    # and a separate feature. `TLslices/atten.env` takes its place: it is a genuine
                    # trapped-mode case with loss in *both* the water and the half-space, and it
                    # exercises the dB/(km·Hz) units that no other case here does.
                    #
                    # Measured 2026-08-09:
                    #   TLslices/atten.env      44 modes  Re 7.2e-7  Im 1.8e-3
                    #   SedAtten/calibK.env     11 modes  Re 1.8e-4  Im 1.1e-2
                    #   SedAtten/calibS_0.6dB   11 modes  Re 2.3e-4  Im 1.5e-2
                    #   SedAtten/calibS_noloss  11 modes  Re 6.0e-10 Im 0 (exactly, both sides)
                    oalib_atten_cases = [
                        (file="TLslices/atten.env", freq=10.0, nmodes=44, kr_rtol=1e-5, α_rtol=1e-2, lossy=true),
                        (file="SedAtten/calibK.env", freq=250.0, nmodes=11, kr_rtol=1e-3, α_rtol=5e-2, lossy=true),
                        (
                            file="SedAtten/calibS_0.6dB.env",
                            freq=250.0,
                            nmodes=11,
                            kr_rtol=1e-3,
                            α_rtol=8e-2,
                            lossy=true,
                        ),
                        (
                            file="SedAtten/calibS_noloss.env",
                            freq=250.0,
                            nmodes=11,
                            kr_rtol=1e-6,
                            α_rtol=1e-12,
                            lossy=false,
                        ),
                    ]

                    @testset "$(case.file)" for case in oalib_atten_cases
                        path = joinpath(oalib_tree, case.file)
                        if !isfile(path)
                            @info "Not present in this Acoustics Toolbox checkout — skipped." path
                            continue
                        end
                        parsed = KR.read_env_file(path)
                        @test is_lossy(parsed.env) == case.lossy

                        c = KR.compare_with_fortran(parsed.env, case.freq)
                        @test c.n_julia == case.nmodes
                        @test c.n_fortran == case.nmodes
                        @test KR.max_kr_reldiff(c) < case.kr_rtol
                        @test KR.min_mode_corr(c) > 0.999
                        @test KR.max_alpha_reldiff(c) < case.α_rtol

                        if case.lossy
                            @test all(c.alpha_julia .< 0)
                            @test all(c.alpha_fortran .< 0)
                        else
                            @test all(iszero, c.alpha_julia)
                            @test all(iszero, c.alpha_fortran)
                        end
                    end

                    # The lossless sibling of calibS is the control for the pair: same waveguide,
                    # same 11 modes, attenuation the only difference. It must agree far more tightly
                    # than the lossy one, and it does -- by six orders of magnitude.
                    lossless = joinpath(oalib_tree, "SedAtten/calibS_noloss.env")
                    lossy = joinpath(oalib_tree, "SedAtten/calibS_0.6dB.env")
                    if isfile(lossless) && isfile(lossy)
                        cl = KR.compare_with_fortran(KR.read_env_file(lossless).env, 250.0)
                        cy = KR.compare_with_fortran(KR.read_env_file(lossy).env, 250.0)
                        @test KR.max_kr_reldiff(cl) < KR.max_kr_reldiff(cy)
                        # ...and the real parts still agree well, because attenuation shifts them
                        # only at second order.
                        @test cl.kr_julia ≈ cy.kr_julia rtol = 1e-3
                    end
                end

                @testset "M5.1: the toolbox's own attenuation cases parse" begin
                    # These are the files the milestone is aimed at, read in place (GPL-3 tree, MIT
                    # package). Each row is what the .env actually declares -- checked by eye against
                    # the file, not against the reader.
                    atten_cases = [
                        (file="SedAtten/calibK.env", units=:dB_per_wavelength, αb=0.5, αp=0.0, freq=250.0),
                        (file="SedAtten/calibS_0.6dB.env", units=:dB_per_wavelength, αb=0.6, αp=0.0, freq=250.0),
                        (file="SedAtten/calibS_noloss.env", units=:dB_per_wavelength, αb=0.0, αp=0.0, freq=250.0),
                        (file="TLslices/atten.env", units=:dB_per_kmHz, αb=0.001, αp=0.001, freq=10.0),
                    ]
                    @testset "$(case.file)" for case in atten_cases
                        path = joinpath(oalib_tree, case.file)
                        if !isfile(path)
                            @info "Not present in this Acoustics Toolbox checkout — skipped." path
                            continue
                        end
                        parsed = KR.read_env_file(path)
                        @test parsed.atten_units === case.units
                        @test parsed.env.αb ≈ case.αb
                        @test all(≈(case.αp), parsed.env.α.α)
                        @test parsed.freqs[1] == case.freq
                        @test is_lossy(parsed.env) == (case.αb > 0 || case.αp > 0)
                    end

                    # The power-law variant is named rather than mis-read as plain dB/m.
                    powlaw = joinpath(oalib_tree, "SedAtten/calibS_PowLaw.env")
                    if isfile(powlaw)
                        err = try
                            KR.read_env_file(powlaw)
                            nothing
                        catch e
                            e
                        end
                        @test err isa KR.UnsupportedEnvFeature
                        @test err.feature == "power-law attenuation"
                    end
                end

                @testset "M6.5: the toolbox's newly readable options against kraken.exe" begin
                    # Every case here was refused by the reader before Milestone 6. CLOW/CHIGH come
                    # from the file, because several decks narrow the band on purpose: MunkK1525 and
                    # gulf_rd stop at 1525 m/s, so Fortran reports only the slow modes while
                    # Kraken.jl finds every trapped one, and the leading modes are what is compared.
                    #
                    # The Munk decks' worst modes are the last trapped ones, 99-102 at the half-space
                    # cutoff. That is not the interpolator: forcing C-linear on both sides gives the
                    # same 6.1e-5 and 0.9944. wedge.env's worst are its grazing modes near kᵣ = 0
                    # (mode 59, kᵣ = 0.0147 at ω/c = 0.105), where relative error is inflated.
                    toolbox_options = [
                        (file="TLslices/pekeris.env", nmodes=44, kr_rtol=1e-6, corr=0.9999),
                        (file="Noise/Pekeris/pekeris.env", nmodes=26, kr_rtol=1e-4, corr=0.9999),
                        (file="MunkLeaky/MunkK1525.env", nmodes=28, kr_rtol=1e-5, corr=0.9999),
                        (file="Gulf/gulf_rd.env", nmodes=63, kr_rtol=1e-5, corr=0.9999),
                        (file="Munk/MunkK.env", nmodes=102, kr_rtol=2e-4, corr=0.99),
                        (file="Munk/MunkS.env", nmodes=102, kr_rtol=2e-4, corr=0.99),
                        (file="wedge/wedge.env", nmodes=59, kr_rtol=1e-2, corr=0.999),
                    ]
                    @testset "$(case.file)" for case in toolbox_options
                        path = joinpath(oalib_tree, case.file)
                        if !isfile(path)
                            @info "Not present in this Acoustics Toolbox checkout — skipped." path
                            continue
                        end
                        parsed = KR.read_env_file(path)
                        c = KR.compare_with_fortran(parsed.env, parsed.freqs[1]; clow=parsed.clow, chigh=parsed.chigh)
                        @test c.n_fortran == case.nmodes
                        @test c.n_julia >= c.n_fortran
                        @test KR.max_kr_reldiff(c) < case.kr_rtol
                        @test KR.min_mode_corr(c) > case.corr
                    end
                    # What each case exercises, so a change in how the reader maps them is caught here.
                    modes_of(file) =
                        let e = KR.read_env_file(joinpath(oalib_tree, file)).env
                            (e.c.mode, e.bottom_bc)
                        end
                    if isfile(joinpath(oalib_tree, "wedge/wedge.env"))
                        @test modes_of("Munk/MunkK.env") == (:n2_linear, AcousticHalfspace())
                        @test modes_of("Munk/MunkS.env") == (:cubic_spline, AcousticHalfspace())
                        @test modes_of("wedge/wedge.env") == (:c_linear, PressureRelease())
                    end
                end

                @testset "M6.4: every file that reads writes back to the same problem" begin
                    # read_env_file -> write_env_file must preserve what the file *means*, and the only
                    # judge of that is kraken.exe: run it on a copy of the original and on the rewrite,
                    # and the two must agree on the mode count and on every wavenumber. The rewrite
                    # keeps the file's NMESH, CLOW/CHIGH and title; the source and receiver depths
                    # only tabulate the modes and are left at the writer's defaults.
                    #
                    # Decks kraken.exe cannot run *as shipped* have no reference to compare with and
                    # are counted, not failed. They are BELLHOP inputs whose NMESH is far too coarse
                    # for KRAKEN ("Mesh is too coarse"), or whose trailing records stop it before it
                    # writes a .mod.
                    function run_deck(path, freq)
                        dir = mktempdir()
                        try
                            cp(path, joinpath(dir, "case.env"))
                            cmd = Cmd(`$(KR.kraken_cmd()) case`; dir=dir)
                            run(pipeline(ignorestatus(cmd); stdout=devnull, stderr=devnull))
                            prt = joinpath(dir, "case.prt")
                            report = isfile(prt) ? read(prt, String) : ""
                            isempty(KR._error_lines(report)) && isfile(joinpath(dir, "case.mod")) || return nothing
                            modes = KR.read_mod_file(joinpath(dir, "case.mod"); freq=freq)
                            grp = try
                                KR.read_grp(prt; freq=freq)
                            catch
                                nothing
                            end
                            return (; kr=KR._best_fortran_kr((; kᵣ=modes.kᵣ, grp)), α=imag.(modes.kᵣ))
                        finally
                            rm(dir; recursive=true, force=true)
                        end
                    end

                    # Two decks take minutes in kraken.exe alone (measured 2026-09-16: 72 s and 167 s,
                    # 626 and 1221 modes) and round-tripped exactly then. `Dickins/Precalc/DickinsK.env`
                    # is the same rigid-bottom waveguide at 1221 modes in 1.6 s and stays in.
                    slow_roundtrips = ("Dickins/March/DickinsK_rd.env", "Dickins/Precalc/DickinsK_rd.env")

                    tally = Dict(:compared => 0, :original_fails => 0)
                    for rel in sort(KR.categorize_env_tree(oalib_tree).supported)
                        rel in slow_roundtrips && continue
                        path = joinpath(oalib_tree, rel)
                        parsed = KR.read_env_file(path)
                        freq = parsed.freqs[1]
                        original = run_deck(path, freq)
                        if original === nothing
                            tally[:original_fails] += 1
                            continue
                        end
                        rewritten = mktempdir() do dir
                            kw = (; title=parsed.title, nmesh=parsed.nmesh, rmax=something(parsed.rmax, 10.0))
                            parsed.clow === nothing || (kw = (; kw..., clow=parsed.clow, chigh=parsed.chigh))
                            run_deck(KR.write_env_file(joinpath(dir, "rt"), parsed.env, freq; kw...), freq)
                        end
                        tally[:compared] += 1
                        ok = rewritten !== nothing && rewritten.kr == original.kr && rewritten.α == original.α
                        ok || @info "Round trip changed the problem" rel parsed.topopt parsed.botopt
                        @test ok
                    end
                    @info "OALIB round trip" compared = tally[:compared] original_fails = tally[:original_fails]
                    # 94 compared, all exact, and 115 unrunnable as shipped (measured 2026-09-16).
                    @test tally[:compared] > 80
                end

                @testset "categorized report over the whole tree" begin
                    report = KR.categorize_env_tree(oalib_tree)
                    @test report.total > 100
                    @test length(report.supported) > 20
                    # Every rejection carries a reason; that list is the Milestone 5/6 backlog.
                    @test !isempty(report.unsupported)
                    @test all(!isempty, values(report.unsupported))
                    # Plain compressional attenuation stopped being a blocker in Milestone 5.1 --
                    # what remains under that name is the power law and the added volume-attenuation
                    # laws (Thorp, Francois-Garrison, biological), which are separate features.
                    @test !haskey(report.unsupported, "attenuation")
                    text = sprint(KR.print_env_tree_report, report)
                    @test occursin("Scanned $(report.total) .env files", text)
                    @test occursin("boundary", text)
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
