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
# Run a single TestItems file
julia --project=test -e 'using TestItemRunner; @run_package_tests filter=t->t.filename=="test/environment_tests.jl"'

# Run a single @testset file
julia --project=test -e 'using Kraken; include("test/numerical_methods_tests.jl")'
```

## Test Dependencies

Required packages (in test/Project.toml):
- `Test`, `TestItems`, `TestItemRunner` - Test frameworks
- `ForwardDiff`, `FiniteDiff` - AD testing
- `BenchmarkTools` - Performance tests
- `Roots` - Root finding methods
