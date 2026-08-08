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

Julia 1.12.6, macOS arm64 (M1), 1 thread inside `Pkg.test()`. Every run below is green: 0 failures,
0 errors, 0 broken.

| Run | End of M1 (`df99dfb`) | End of M2 (`ce38000`) |
|---|---|---|
| `Pkg.test()` | **258** / 2m05s | **380** / 2m05s |
| `KRAKEN_RUN_PERFORMANCE_TESTS=true Pkg.test()` | **282** / 2m25s | **404** / 2m20s |

Per file, so a silent drop in coverage is visible in a diff:

| File | End of M1 | End of M2 |
|---|---|---|
| `environment_tests.jl` | 39 | 161 |
| `integration_tests.jl` | 98 | 98 |
| `numerical_methods_tests.jl` | 73 | 73 |
| `automatic_differentiation_tests.jl` | 48 | 48 |
| `performance_tests.jl` (opt-in) | 24 | 24 |

The M2 jump is the B1–B5 regression tests added in task 2.5; the B4 bisection sweep (8 environments
× 4 frequencies) is 82 of the 122 new assertions on its own.

The suite also runs on Julia 1.10, the declared compat lower bound — verified directly, not assumed
(see the note in the plan's Architecture Decisions).

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
| `fortran_reference_tests.jl` | Test | Cross-validation against unmodified `kraken.exe`; self-skips when no binary resolves |

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

`fortran_reference_tests.jl` runs the suite against **unmodified** Fortran KRAKEN. The harness lives
in `test/reference/KrakenReference.jl` — test-only, deliberately not part of the package, because
validation machinery is not public API.

Binaries come from `AcousticsToolbox_jll`, which ships prebuilt `kraken.exe`/`krakenc.exe` for every
platform, so **CI needs no Fortran toolchain**. To compare against a local Acoustics Toolbox build
instead, point `KRAKEN_FORTRAN_BIN` at a directory containing `kraken.exe`:

```bash
KRAKEN_FORTRAN_BIN=~/programs/AcousticsToolboxOALIB/bin julia --project=. -e 'using Pkg; Pkg.test()'
```

The override wins when it resolves and falls back to the jll when it does not, so a stale value
cannot break the suite. Everything is gated on `KrakenReference.fortran_available()`: a platform with
neither source skips with a message rather than erroring, which is why the file is included
unconditionally from `runtests.jl`.

> **`kraken.exe` exits 0 even on a fatal error.** It writes `STOP Fatal Error: Check the print file
> for details` to stderr and stops with a zero status. Anything driving it must scan the generated
> `.prt` file for `ERROR`, never trust the exit code.

The old `fortran_interface_tests.jl` and its `KRAKEN_RUN_FORTRAN_TESTS` switch were removed in plan
task 1.4: they called an `EnvKRAKEN` API that exists in no module, so they could never have run.

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
