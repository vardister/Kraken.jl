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

    if KR.fortran_available()
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
