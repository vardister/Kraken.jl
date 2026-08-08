# Kraken.jl Test Suite

## Quick Start

Run from the **repo root**, against the **root** environment (`--project=.`). Pkg picks up
`test/Project.toml` automatically as the test environment.

```bash
# Run all tests
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'

# Run with performance tests
KRAKEN_RUN_PERFORMANCE_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```

## Baseline

Recorded at the end of Milestone 1 (commit `df99dfb`, 2026-08-08) on Julia 1.12.6, macOS
arm64 (M1), 1 thread inside `Pkg.test()`. Both runs are green: 0 failures, 0 errors, 0 broken.

| Run | Tests | Suite time |
|---|---|---|
| `Pkg.test()` | **258** | 2m05s |
| `KRAKEN_RUN_PERFORMANCE_TESTS=true Pkg.test()` | **282** | 2m25s |

Per file, so a silent drop in coverage is visible in a diff:

| File | Tests |
|---|---|
| `environment_tests.jl` | 39 |
| `integration_tests.jl` | 98 |
| `numerical_methods_tests.jl` | 73 |
| `automatic_differentiation_tests.jl` | 48 |
| `performance_tests.jl` (opt-in) | 24 |

The suite time is dominated by `automatic_differentiation_tests.jl` (~70 s: it runs `kraken_jl`
under `ForwardDiff` *and* `FiniteDiff` across dozens of parameter points). Add roughly 4 minutes
of dependency precompilation on the first run after a `Manifest.toml` change.

> **Do not use `--project=test` with `Pkg.test()`.** Pkg would then treat the test environment as the
> package under test and fail with `The Project.toml of the package being tested must have a name and
> a UUID entry`. `test/Project.toml` has no `name`/`uuid` on purpose — that is what makes it a test
> environment rather than a package. Correspondingly, the root `Project.toml` has no
> `[extras]`/`[targets]`: the two testing mechanisms are mutually exclusive and `test/Project.toml` wins.

## Test Files

### Core Tests (always run)

| File | Framework | Description |
|------|-----------|-------------|
| `environment_tests.jl` | TestItems | Environment creation, sound speed/density profiles |
| `integration_tests.jl` | TestItems | Full `kraken_jl` workflow, physics validation |
| `numerical_methods_tests.jl` | Test | Core algorithms: `det_sturm`, `bisection`, `find_kr`, `inverse_iteration` |
| `automatic_differentiation_tests.jl` | Test | ForwardDiff compatibility, group speed calculations |

### Optional Tests

| File | Env Variable | Description |
|------|--------------|-------------|
| `performance_tests.jl` | `KRAKEN_RUN_PERFORMANCE_TESTS=true` | Benchmarks, memory usage, scaling |

Every timing in `performance_tests.jl` is measured after a warm-up call, so compilation is excluded.
Thresholds are set several times above the measured value so they catch order-of-magnitude
regressions without flaking on slower CI runners. The frequency- and depth-scaling ratios are
*reported*, not asserted — wall-clock ratios between sub-second runs are dominated by timer noise —
and those testsets assert a generous absolute ceiling plus the mode counts instead.

### Script Files (manual execution)

| File | Description |
|------|-------------|
| `timings_vs_fortran.jl` | Manual timing comparisons (requires DrWatson; calls the removed `EnvKRAKEN` API, so it does not currently run) |

### Fortran cross-validation

There is none at present. `fortran_interface_tests.jl` and the `KRAKEN_RUN_FORTRAN_TESTS` switch were
removed in plan task 1.4: they called an `EnvKRAKEN` API that exists in no module, so they could never
have run. Milestone 3 replaces them with a `test/reference/` harness that drives unmodified
`kraken.exe` from `AcousticsToolbox_jll` over `.env`/`.mod` files.

Prior Enzyme.jl AD experiments lived in `test/ad_tests.jl`, deleted in the same task. They are
recoverable if Milestone 4 wants them: `git show 580649c:test/ad_tests.jl`.

## Running Individual Test Files

These run *inside* the test environment (`--project=test`), which is fine — the restriction above
applies only to `Pkg.test()`. They do require a one-time setup step:

```bash
# ONE-TIME: dev-link the test environment to the local package
julia --project=test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'

# Verify it points at the working tree, not ~/.julia/packages/Kraken/...
julia --project=test -e 'using Kraken; println(pathof(Kraken))'
```

> **Why this is required, and not optional.** `Kraken` is registered in the General registry, and
> `test/Manifest.toml` is gitignored. On a fresh clone, `--project=test` therefore resolves `Kraken`
> from the *registry* — you would be testing the last released version instead of your working tree,
> silently and with no error. The `Pkg.develop` call above overrides that with a path entry.
> `Pkg.test()` from the root environment is immune to this: it always uses the local package.

```bash
# Run a single TestItems file. `t.filename` is an ABSOLUTE path, so match with `endswith`
# — a `==` against "test/environment_tests.jl" silently selects zero test items and reports
# a green "Package | 0 total".
julia --project=test -e 'using TestItemRunner; @run_package_tests filter=t->endswith(t.filename, "environment_tests.jl")'

# Run a single @testset file
julia --project=test -e 'using Kraken; include("test/numerical_methods_tests.jl")'
```

## Test Dependencies

Required packages (in test/Project.toml):
- `Test`, `TestItems`, `TestItemRunner` - Test frameworks
- `ForwardDiff`, `FiniteDiff` - AD testing
- `BenchmarkTools` - Performance tests
- `Roots` - Root finding methods
