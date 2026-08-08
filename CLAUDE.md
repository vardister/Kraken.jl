# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Kraken.jl is a pure-Julia reimplementation of Michael Porter's KRAKEN normal-mode model for underwater
acoustic propagation (from the Acoustics Toolbox), plus ideas from UnderwaterAcoustics.jl. It computes
horizontal wavenumbers and mode shapes for a range-independent underwater environment via a
finite-difference discretization of the depth-separated wave equation, then uses those modes to
synthesize the acoustic pressure field. It is cross-validated against unmodified Fortran KRAKEN by a
test-only harness under `test/reference/` (see below). There is no Fortran source, shared library, or
`Makefile` in this repo, and nothing here links against Fortran code.

## Commands

> **Never run Julia through Bash. Use the `kaimon` MCP.** This is not a preference about speed — the
> user watches the shared REPL, and a `julia -e ...` subprocess is invisible to them. It applies to
> *every* Julia invocation, with no exceptions, including the full test suite:
>
> | Instead of | Use |
> |---|---|
> | `julia --project=. -e 'using Pkg; Pkg.test()'` | `run_tests(project_path="<repo root or worktree>")` |
> | `julia -e '...'`, `julia --project=... -e '...'` | `ex(e="...", ses="<8-char key>")` |
> | `julia -e 'using JuliaFormatter; format(...)'` | `format_code(path="...")` |
> | `julia --project=... -e 'Pkg.add(...)'` | `pkg_add(packages=[...])` |
>
> No session connected for this checkout? Start one — `start_session(project_path="…")` — rather than
> falling back to Bash. In a git worktree, point it at the *worktree* path, or Revise reloads the
> wrong copy of `src/`. Env-var-prefixed runs (`KRAKEN_RUN_PERFORMANCE_TESTS=true …`) go through the
> MCP too: set them with `ex(e="ENV[\"KRAKEN_RUN_PERFORMANCE_TESTS\"] = \"true\"")` in the session
> that will run the suite. Bash stays for git, file moves, and non-Julia tooling only.
>
> The `julia ...` command lines below are the *documented invocations* — what CI runs and what a
> human types in a terminal. They are the reference for which project flag is correct; they are not
> an instruction to shell out.

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

# Cross-validate against a local Acoustics Toolbox build instead of the AcousticsToolbox_jll
# binaries, and enable the test cases built from the toolbox's own .env files. Both optional —
# without them the suite uses the jll and skips the toolbox cases (which is the CI path).
KRAKEN_FORTRAN_BIN=~/programs/AcousticsToolboxOALIB/bin \
KRAKEN_OALIB_TESTS=~/programs/AcousticsToolboxOALIB/tests \
  julia --project=. -e 'using Pkg; Pkg.test()'

# ONE-TIME setup, REQUIRED before the single-file invocations below.
# Run it once per clone AND once per git worktree — `.gitignore`'s bare `Manifest.toml` pattern
# matches at every depth, so test/Manifest.toml is never carried into a new worktree either.
# `Kraken` is registered in General, so without this dev-link `--project=test` silently resolves
# Kraken to the RELEASED version from the registry instead of your working tree, and you test code
# you did not write. (`Pkg.test()` from the root env is immune — it always uses the local package,
# which is why the two run modes can disagree.) `test/runtests.jl` now checks this on startup and
# fails with the fix command rather than letting the suite run against the wrong package.
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Verify the link points at the working tree, not ~/.julia/packages/Kraken/...
julia --project=test -e 'using Kraken; println(pathof(Kraken))'

# Run a single TestItems-based file (environment_tests.jl, integration_tests.jl).
# `t.filename` is an ABSOLUTE path — match with `endswith`, not `==`, or the filter selects
# zero test items and reports a green "Package | 0 total".
julia --project=test -e 'using TestItemRunner; @run_package_tests filter=t->endswith(t.filename, "environment_tests.jl")'

# Run a single @testset-based file (numerical_methods_tests.jl, automatic_differentiation_tests.jl)
julia --project=test -e 'using Kraken; include("test/numerical_methods_tests.jl")'

# Format (JuliaFormatter, config in .JuliaFormatter.toml: indent=4, style=blue, margin=120).
# Formats in place AND reports — this is both the fix command and exactly what CI runs, so run it
# before pushing. Do NOT use `format(".")`: JuliaFormatter ignores .gitignore and descends into
# .git/, where this repo has branch refs whose names end in `.jl`.
julia .github/format_check.jl

# Build the docs. CI builds them too, so this must stay green.
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'   # one-time
julia --project=docs docs/make.jl
```

**Julia 1.10 is the declared lower bound and it is real** — CI runs the suite on it. Two things bite:
`@views x[1:(end-1)] .-= y` (an updating broadcast under `@views`) is a *syntax error* on 1.10, so
write `@views x[1:(end-1)] .= x[1:(end-1)] .- y` instead. Verify locally with
`juliaup add 1.10 && julia +1.10 --startup-file=no --project=. -e 'using Pkg; Pkg.test()'` — note
`--startup-file=no`, or the user's `startup.jl` fails on a missing Revise in that depot.

Test files fall into three categories (see `test/README.md` for the up-to-date table):
- **TestItems** files (`environment_tests.jl`, `integration_tests.jl`) — use `@testitem`, run via `@run_package_tests` in `runtests.jl`.
- **Test** files (`numerical_methods_tests.jl`, `automatic_differentiation_tests.jl`, `fortran_reference_tests.jl`, `performance_tests.jl`) — plain `@testset`, `include`d from `runtests.jl`.
- **Manual scripts** (`timings_vs_fortran.jl` needs DrWatson) — not wired into `runtests.jl`, run by hand.

**Anything added to `test/` must be declared in `test/Project.toml`, stdlibs included.**
`--project=test` resolves against a manifest that already contains every stdlib indirectly, so
`using LinearAlgebra` succeeds there even when it was never declared; `Pkg.test()` builds from the
declared dependencies only and fails. A green single-file run proves less than it looks like it does
— run the full `Pkg.test()` before pushing.

### Fortran cross-validation (`test/reference/`)

`test/reference/KrakenReference.jl` is a test-only module — validation machinery is not public API,
so it is deliberately not in `src/`. It `include`s, in order:

| File | Provides |
|---|---|
| `env_writer.jl` | `write_env_file` / `env_file_string` — an `UnderwaterEnv` → KRAKEN `.env` |
| `mod_reader.jl` | `read_mod_file` (binary `.mod`), `read_grp` (the `.prt` Group Speed table), `has_group_speeds` |
| `env_reader.jl` | `read_env_file` — `.env` → `UnderwaterEnv`; `categorize_env_tree` for bulk feature triage |
| `runner.jl` | `run_fortran_kraken` — write, invoke, read back; raises `FortranKrakenError` |
| `compare.jl` | `compare_with_fortran` — runs both solvers, returns a printable `FortranComparison` |

Binaries come from `AcousticsToolbox_jll` (no Fortran toolchain needed, including in CI). Two env
vars: `KRAKEN_FORTRAN_BIN` points at a directory containing `kraken.exe` to use a local Acoustics
Toolbox build instead; `KRAKEN_OALIB_TESTS` points at a toolbox `tests/` tree to enable the cases
built from its `.env` files (those files are GPL-3 and this package is MIT, so they are read in place,
never vendored). Both are optional — everything is gated on `KrakenReference.fortran_available()` and
skips with a message rather than erroring.

Things that will bite you here, all established by running the code and recorded in
`test/README.md` and the plan's Architecture Decisions:

- **`kraken.exe` exits 0 even on a fatal error.** Scan the generated `.prt` for `ERROR`; never trust
  the exit code.
- **`AcousticsToolbox_jll` v2025.9's `kraken.exe` reports every group speed as `0.00000`.** Its
  wavenumbers and mode shapes are correct; only `VG` is lost, and `krakenc.exe` is unaffected. Use
  `complex=true` for anything needing group speeds.
- **Neither output file dominates on precision.** The `.mod` is single precision throughout; the
  `.prt` prints `Re(kᵣ)` with 10 digits but `alpha` with 2. Wavenumber comparisons want the `.prt`,
  attenuation work wants the `.mod`.
- **The `.prt` group-speed table is subsampled** (`DO mode = 1, M, MAX(1, M/30)`), so `read_grp`
  returns mode indices, not just values.
- **Generated `.env` files use `'C'` (C-linear), not the `'S'` of the checked-in samples**, because
  `SampledSSP` interpolates linearly — `'S'` would have the Fortran solve a different problem.
- **`.env` densities are g/cm³; Kraken.jl's are kg/m³.** The writer and reader convert.

The old `fortran_interface_tests.jl` and its `KRAKEN_RUN_FORTRAN_TESTS` switch were deleted in plan
task 1.4 (they called an `EnvKRAKEN` API that exists nowhere) and are fully replaced by the above.

**KrakenFortran.jl is a separate, optional package and Kraken.jl does not depend on it.** It offers
an *in-process* `ccall` path to Fortran KRAKEN, which is genuinely faster for broadband sweeps. It is
deliberately not the correctness oracle: its sources are a MEX-adapted fork of an older KRAKEN, so
validating against it would prove nothing about agreement with upstream. Do not add it as a
dependency or reference it from the test harness.

## Architecture

`src/Kraken.jl` is the module entry point and `include`s three files, in dependency order:

```julia
include("kraken_core.jl")                  # main FD solver — most work happens here
include("kraken_pekeris.jl")               # closed-form Pekeris (2-layer) analytic model
include("kraken_standard_environments.jl") # canned test environments (pekeris_env, one_layer_env, munk_env, ...)
```

`src/` now contains exactly the files `Kraken.jl` includes, and nothing else. Code that is not part
of the package lives outside it: **`dev/kraken_broadband.jl`** is a standalone script you
`include(...)` by hand (broadband pulse synthesis from mode sums). It is in `dev/` rather than
`src/` so its "not part of the package" status is structural instead of a comment — Milestone 7
promotes it properly, which means moving it back into `src/`, adding the `include` to `Kraken.jl`,
and adding its deps to `Project.toml`.

Plotting lives in **`ext/KrakenMakieExt.jl`**, a package extension triggered by `Makie` (a
`[weakdeps]` entry, so it costs nothing unless you ask for it). `using CairoMakie` or `using GLMakie`
both pull in Makie and activate it. `plot_modes` and `plot_ssp` are exported stubs in `src/Kraken.jl`;
without a backend loaded they raise an actionable error rather than a `MethodError`. The old
`src/kraken_plots.jl` is gone — it took the `Dict` returned by the removed Fortran `kraken` function,
not a `NormalModeSolution`, so it could not have run against the current API.

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
