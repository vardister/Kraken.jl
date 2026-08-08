using Test
using Kraken

include("reference/KrakenReference.jl")
const KR = KrakenReference

# Cross-validation against unmodified Fortran KRAKEN. Plan tasks 3.2-3.7 extend this file with the
# .env writer, the .mod/.prt readers, the runner, and the per-environment comparison suite; for now
# it covers binary resolution only.
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
end
