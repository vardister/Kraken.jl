# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Kraken.jl is a pure-Julia reimplementation of Michael Porter's KRAKEN normal-mode model for underwater
acoustic propagation (from the Acoustics Toolbox), plus ideas from UnderwaterAcoustics.jl. It computes
horizontal wavenumbers and mode shapes for a range-independent underwater environment via a
finite-difference discretization of the depth-separated wave equation, then uses those modes to
synthesize the acoustic pressure field. A prebuilt Fortran shared library (`src/kraken.dylib` /
`.so` / `.dll`, built from `src_fortran/` via `Makefile`/`make`) is used only for cross-validation
against the original Fortran KRAKEN — that comparison layer lives in a separate package,
`KrakenFortran.jl`, not in this repo.

## Commands

Run all commands from the repo root. The package env is `Project.toml` at the root; the test env is
`test/Project.toml`, which Pkg picks up **automatically** when `Pkg.test()` is run from the *root*
environment. `test/Project.toml` deliberately has no `name`/`uuid` (that is correct for a test
environment) and the root `Project.toml` has no `[extras]`/`[targets]` — the two mechanisms are
mutually exclusive and `test/Project.toml` wins. Do **not** run `Pkg.test()` with `--project=test`:
Pkg then treats the test env as the package under test and errors with
`The Project.toml of the package being tested must have a name and a UUID entry`.

```bash
# Instantiate + run the full test suite  (note: --project=. , not --project=test)
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Include performance benchmarks (skipped by default, can be slow)
KRAKEN_RUN_PERFORMANCE_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'

# ONE-TIME setup, REQUIRED before the single-file invocations below.
# `Kraken` is registered in General, and test/Manifest.toml is gitignored — so on a fresh
# clone `--project=test` silently resolves Kraken to the RELEASED version from the registry
# instead of your working tree, and you test code you did not write. This dev-links it.
# (`Pkg.test()` from the root env is immune — it always uses the local package.)
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Verify the link points at the working tree, not ~/.julia/packages/Kraken/...
julia --project=test -e 'using Kraken; println(pathof(Kraken))'

# Run a single TestItems-based file (environment_tests.jl, integration_tests.jl)
julia --project=test -e 'using TestItemRunner; @run_package_tests filter=t->t.filename=="test/environment_tests.jl"'

# Run a single @testset-based file (numerical_methods_tests.jl, automatic_differentiation_tests.jl)
julia --project=test -e 'using Kraken; include("test/numerical_methods_tests.jl")'

# Format (JuliaFormatter, config in .JuliaFormatter.toml: indent=4, style=blue, margin=120)
julia -e 'using JuliaFormatter; format(".")'

# Build the Fortran shared library from src_fortran/ (only needed for Fortran-comparison work)
make
```

Test files fall into three categories (see `test/README.md` for the up-to-date table):
- **TestItems** files (`environment_tests.jl`, `integration_tests.jl`) — use `@testitem`, run via `@run_package_tests` in `runtests.jl`.
- **Test** files (`numerical_methods_tests.jl`, `automatic_differentiation_tests.jl`, `performance_tests.jl`) — plain `@testset`, `include`d from `runtests.jl`.
- **Manual scripts** (`timings_vs_fortran.jl` needs DrWatson) — not wired into `runtests.jl`, run by hand.

There is deliberately no Fortran comparison layer right now. The old `fortran_interface_tests.jl`
and the `KRAKEN_RUN_FORTRAN_TESTS` switch were deleted in plan task 1.4 because they called an
`EnvKRAKEN` API that no longer exists anywhere. Milestone 3 replaces them with `test/reference/`,
driving unmodified `kraken.exe` from `AcousticsToolbox_jll` over `.env`/`.mod` files.

## Architecture

`src/Kraken.jl` is the module entry point and only `include`s three files, in dependency order:

```julia
include("kraken_core.jl")                 # main FD solver — most work happens here
include("kraken_basic.jl")                # currently empty
include("kraken_pekeris.jl")               # closed-form Pekeris (2-layer) analytic model
include("kraken_standard_environments.jl") # canned test environments (pekeris_env, one_layer_env, munk_env, ...)
```

`src/kraken_broadband.jl` and `src/kraken_plots.jl` exist in the tree but are **not** `include`d from
`Kraken.jl` — they're standalone scripts you `include(...)` manually when needed (broadband pulse
synthesis from mode sums, and CairoMakie mode plots respectively). If you promote either to a real
part of the package, add the `include` to `Kraken.jl` and add its deps to `Project.toml`.

### Solve pipeline (`kraken_core.jl`)

1. **Environment**: `UnderwaterEnv` wraps a `SampledSSP` (sound-speed profile) and `SampledDensity`
   profile plus bottom half-space properties (`cb`, `ρb`) and layer depths. Built from KRAKEN-file-style
   `ssp`/`layers`/`sspHS` matrices (see `kraken_standard_environments.jl` for the matrix layout) via
   `UnderwaterEnv(ssp, layers, sspHS)`.
2. **Discretization**: `AcousticProblemProperties(env, freq)` picks a per-layer mesh (`Nz_vec`, `Δz_vec`,
   `zn_vec`) sized to ~20 points per wavelength at the fastest bottom sound speed.
3. **Cache**: `AcousticProblemCache(env, props)` builds the finite-difference tridiagonal system
   (`a_vec`/`e_vec`/`A::Tridiagonal`) that root-finding and inverse iteration reuse and mutate in place.
4. **Root finding**: `det_sturm` evaluates the Sturm sequence (mode-counting determinant) for a trial
   wavenumber `kr`. `bisection` uses Sturm sequence mode counts to bracket each mode's wavenumber, then
   `find_kr`/`solve_for_kr` refine each bracket via NonlinearSolve.jl's `ITP()` (bracketing solver on
   `IntervalNonlinearProblem`).
5. **Mode shapes**: `inverse_iteration` recovers the eigenvector (mode shape) for each converged `kr`
   by inverse-iterating the cached tridiagonal system, then normalizes it against the water column +
   bottom half-space energy.
6. **Convergence**: `kraken_jl(env, freq)` is the top-level entry point — it repeats steps 2-5 across
   successively finer meshes and Richardson-extrapolates (`richard_extrap`) the wavenumbers-squared
   across mesh levels until convergence (or `n_meshes` is exhausted), returning a `NormalModeSolution`
   (`kr`, `modes`, `env`, `props`).

`kraken_pekeris.jl` is a separate, independent solve path: for the 2-layer Pekeris waveguide it solves
the transcendental dispersion relation directly (`find_kr(::PekerisUnderwaterEnv, freq)` via bracketing
+ `Roots.jl`) rather than going through the finite-difference/Sturm-sequence machinery — useful as a
ground truth / cross-check for the general FD solver.

### Automatic differentiation

The solver is written to support `ForwardDiff.jl` differentiation straight through (`find_kr` →
`kraken_jl`), e.g. to get group speeds as `dω/dkr` (see README "Calculating group speeds"). Keep new
numerical code AD-friendly (avoid non-differentiable branching on Dual numbers, prefer `NaNMath`-style
functions already used in `kraken_core.jl`/`kraken_pekeris.jl` for `sqrt`/`tan` etc. to avoid
domain errors under AD perturbation).

### Environment construction

`kraken_standard_environments.jl` mirrors the Acoustics-Toolbox `.env`-file structure: `ssp` is an
`[depth, cp, cs, ρ, αp, αs]` matrix, `layers` marks per-layer depth boundaries, and `sspHS` gives the
surface/bottom half-space properties — this is the shape `UnderwaterEnv`/`UnderwaterEnvFORTRAN`
constructors expect. New standard-environment helpers should follow the same matrix shape so they work
with both the Julia and Fortran-interface code paths.
