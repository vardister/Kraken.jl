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

| Run | End of M1 (`df99dfb`) | End of M2 (`ce38000`) | End of M3 |
|---|---|---|---|
| `Pkg.test()` | **258** / 2m05s | **380** / 2m05s | **925** / 2m53s |
| `KRAKEN_RUN_PERFORMANCE_TESTS=true Pkg.test()` | **282** / 2m25s | **404** / 2m20s | **949** / 3m10s |

The M3 numbers are with an Acoustics Toolbox checkout present. Without one — which is the CI case —
the toolbox cases skip and `Pkg.test()` reports **898**, still green.

Per file, so a silent drop in coverage is visible in a diff:

| File | End of M1 | End of M2 | End of M3 |
|---|---|---|---|
| `environment_tests.jl` | 39 | 161 | 161 |
| `integration_tests.jl` | 98 | 98 | 98 |
| `numerical_methods_tests.jl` | 73 | 73 | 73 |
| `automatic_differentiation_tests.jl` | 48 | 48 | 48 |
| `fortran_reference_tests.jl` | — | — | 545 |
| `performance_tests.jl` (opt-in) | 24 | 24 | 24 |

The M2 jump is the B1–B5 regression tests added in task 2.5; the B4 bisection sweep (8 environments
× 4 frequencies) is 82 of the 122 new assertions on its own. The M3 jump is the Fortran
cross-validation harness, which also accounts for the ~50 s of added wall time — two thirds of it in
the `munk` sweep alone, which solves for up to 817 modes at 400 Hz.

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
| `reverse_ad_tests.jl` | Test | The `ChainRulesCore` rules in `src/kraken_ad.jl`: each rule in isolation, then Zygote and Mooncake against ForwardDiff and central differences across three environments and three targets, plus `ChainRulesTestUtils.test_rrule` |
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

#### Measured agreement

The numbers the cross-validation suite asserts against. Measured 2026-08-08 on macOS arm64, five
standard environments at 25 / 50 / 100 / 200 / 400 Hz — **identical against `AcousticsToolbox_jll`
and against a local 2023 OALIB build**, so these are properties of the solver, not of one binary.

`max rel Δkᵣ` is the largest relative wavenumber difference over all modes and all five frequencies;
`min corr` is the smallest mode-shape correlation. Update this table when the tolerances in
`fortran_reference_tests.jl` change, so a loosened bound is visible in a diff.

| Environment | Modes found (25→400 Hz) | max rel Δkᵣ | min corr | tolerance asserted |
|---|---|---|---|---|
| `pekeris_env` | 1, 2, 5, 9, 19 | 2.1e-8 | 0.9999997 | 1e-6 / 0.9999 |
| `one_layer_env` | 1, 3, 5, 10, 21 | 1.4e-7 | 0.9999996 | 1e-6 / 0.9999 |
| `one_layer_slope_env` | 1, 2, 5, 10, 21 | 3.8e-6 | 0.9999994 | 2e-5 / 0.9999 |
| `two_layer_slope_env` | 2, 5, 10, 19, 39 | 2.7e-5 | 0.9999484 | 1e-4 / 0.9995 |
| `munk_env` | 51, 102, 204, 409, 817 | 6.4e-6 | 0.9999720 | 5e-5 / 0.9999 |

Error grows with the number of gradient layers, not with frequency — `two_layer_slope_env` is the
worst case at its *lowest* frequency, where the mesh is coarsest relative to the mode structure.

#### Acoustics Toolbox test cases

`test/reference/env_reader.jl` parses a KRAKEN `.env` file back into an `UnderwaterEnv`, which turns
the environments shipped with the Acoustics Toolbox into test cases. Those files are **GPL-3 while
this package is MIT**, so they are read in place rather than copied into this repo:

```bash
KRAKEN_OALIB_TESTS=~/programs/AcousticsToolboxOALIB/tests julia --project=. -e 'using Pkg; Pkg.test()'
```

Without that tree the toolbox cases skip; the reader itself stays covered by round-tripping this
repo's own `test/standard_envs/` files through `write_env_file` → `read_env_file`.

Coverage as of 2026-08-08 — **65 of 402 `.env` files are supported today, and all 65 cross-validate**
(max relative Δkᵣ 3.6e-6, min correlation 0.999992, mode counts matching within one):

| Blocker | Files | Unblocked by |
|---|---|---|
| attenuation | 114 | Milestone 5 |
| top boundary (not vacuum) | 65 | Milestone 6 |
| bottom boundary (not an acoustic half-space) | 50 | Milestone 6 |
| SSP interpolation over a varying profile | 28 | Milestone 6 |
| added volume attenuation (THORP / Francois-Garrison / biological) | 27 | Milestone 5 |
| bottom half-space is not the fastest medium | 19 | leaky modes (M5 stretch) |
| elastic layer | 7 | out of scope |
| interfacial roughness | 5 | out of scope |
| analytic SSP; profile not starting at the surface | 2 | — |
| not a KRAKEN deck (BELLHOP3D `'H'`/`'Q'` SSP options, malformed) | 20 | n/a |

Regenerate this with `KrakenReference.categorize_env_tree` and `print_env_tree_report`. A file that
uses an unsupported feature is *named*, never approximated — the whole point is that a case Kraken.jl
cannot model fails with "unsupported feature: attenuation" rather than silently mis-parsing into a
plausible environment that then disagrees with Fortran for reasons nobody can find.

Two caveats the suite encodes rather than hides:

- **Mode counts may differ by one at cutoff.** `bisection` searches phase speeds up to `0.9999·cb`
  while the generated `.env` asks KRAKEN for up to `cb` exactly, so Kraken.jl's window is marginally
  narrower and can miss a mode sitting right at the bottom cutoff. Seen on `munk_env` at 100 Hz
  (204 vs 205) and 400 Hz (817 vs 818), always with Julia one short. Only that environment is
  allowed any slack; anywhere else a count mismatch fails.
- **`AcousticsToolbox_jll`'s `kraken.exe` reports every group speed as `0.00000`.** Wavenumbers and
  mode shapes are correct and match the 2023 OALIB build digit for digit; only `VG` is lost, and the
  same jll's `krakenc.exe` is unaffected. `compare_with_fortran(...; group_speeds=true)` therefore
  re-runs with `krakenc.exe` to get a reference. Group speeds are off by default because obtaining
  the Julia side means a ForwardDiff pass through the whole solver (~4 s).

### AD validated against Fortran (plan task 4.7)

`test/reverse_ad_tests.jl` checks the Milestone 4 rules against ForwardDiff and FiniteDiff, but all of
those differentiate the same `det_sturm` — an error in the determinant moves every one of them
together. Two checks in `fortran_reference_tests.jl` use `kraken.exe` as the oracle instead.

**Group speeds** — KRAKEN computes them numerically and prints them; Kraken.jl differentiates for
them. Max relative difference at 100 Hz, on a pinned `nmesh = 4000` (measured 2026-08-09):

| environment | max rel. difference |
|---|---|
| `pekeris` | 1.8e-6 |
| `one_layer` | 2.6e-6 |
| `one_layer_slope` | 3.0e-6 |
| `two_layer_slope` | 1.2e-4 |
| `munk` | 8.8e-5 |

**Pinning the Fortran mesh matters here and nowhere else so far.** On KRAKEN's automatic mesh
(`NMESH = 0`) `two_layer_slope` disagrees by 3.4e-3 — above the 0.1% bound — because the auto mesh is
too coarse across its 20 m layers to give an accurate numerical `dω/dkᵣ`. Tightening Kraken.jl's own
tolerances changes that number in the sixth digit; refining Fortran's mesh drops it 29×. The error was
Fortran's discretization, not Kraken.jl's.

**Gradients against a finite-differenced `kraken.exe`** — the sharper check, since it works for *any*
parameter rather than only frequency. Perturb one parameter, write two `.env` files, run the binary on
each, central-difference `Re(kᵣ)`. Zygote vs that oracle, mode 1 at 100 Hz:

| environment | parameter | step | rel. difference |
|---|---|---|---|
| `pekeris` | `c0` (sound speed) | 1e-3 | 6.1e-7 |
| `pekeris` | `ρ0` (density) | 1e-2 | 6.1e-5 |
| `pekeris` | `depth` (thickness) | 1e-3 | 1.8e-5 |
| `pekeris` | `cb` (control) | 1e-3 | 1.2e-4 |
| `one_layer` | `c1` (sound speed) | 1e-3 | 4.0e-4 |
| `one_layer` | `ρ1` (density) | 1e-3 | 1.8e-4 |
| `one_layer` | `h1` (thickness) | 1e-2 | 2.1e-3 |
| `one_layer` | `c0` (control) | 1e-3 | 2.2e-6 |

The step is per-parameter because the right one is set by the size of the derivative, not the
parameter. The `.prt` prints ten digits of `kᵣ ≈ 0.42`, so a difference below ~1e-10 is quantization:
`∂kᵣ/∂h1` is 3.7e-8, and at `h = 1e-3` the two runs differ by ~14 units in the last printed place —
a 4.6% error that falls to 0.21% at `h = 1e-2`. Stepping larger is not free either: `cb` at `h = 1e-1`
puts the half-space below the water column and the writer rejects the environment.

Reverse mode is separately required to reproduce forward mode, measured **against the gradient's own
scale** rather than entrywise (pekeris 9.6e-12, one_layer 2.1e-11). These gradients span four orders
of magnitude — `∂kᵣ/∂h1` is 1.3e-4 of the largest entry — so an entrywise bound on the smallest
components demands agreement finer than either method achieves. An entrywise `rtol = 1e-8` fails on
`h1` at 1.3e-7 while the two agree to 2.1e-11 on scale. Same reasoning as `relerr_norm` in
`reverse_ad_tests.jl`.

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
> `test/Manifest.toml` is gitignored — `.gitignore` carries a bare `Manifest.toml` pattern, which
> matches at *every* depth, so the rule aimed at the root manifest catches the test one too. On a
> fresh clone — **or in a new git worktree**, which is the easy case to forget, since ignored files
> are not carried over — `--project=test` therefore resolves `Kraken` from the *registry*. You would
> be testing the last released version instead of your working tree, silently and with no error; it
> surfaces as `UndefVarError` on every symbol added since that release. The `Pkg.develop` call above
> overrides that with a path entry. `Pkg.test()` from the root environment is immune: it always uses
> the local package, which is why the two run modes can disagree — and the immune one is not the one
> to trust when they do.
>
> `test/runtests.jl` guards this on startup: `check_testing_this_checkout()` compares
> `pkgdir(Kraken)` against the directory `test/` lives in and aborts with the fix command if they
> differ, so the mistake can no longer be mistaken for a code failure.
>
> **The manifest stays gitignored on purpose.** Committing it would pin the dev-link for every
> clone, but Kraken.jl is a library: an un-pinned test environment is what makes CI re-resolve and
> tell you when a new `DataInterpolations`/`NonlinearSolve` release breaks the package. A committed
> manifest also couples the tree to one Julia version (CI runs a 1.10 leg) and only works while the
> `[[deps.Kraken]]` entry stays relative — `Pkg.develop` from outside the repo writes an *absolute*
> path, which would hard-code one developer's machine into the repo. The startup guard is the fix
> that has none of those costs. (Once the `julia = "1.10"` bound is dropped, a `[sources]` entry in
> `test/Project.toml` — Pkg 1.11+ — would make the dev-link declarative and tracked; it is not
> usable while the LTS is supported.)

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
- `AcousticsToolbox_jll` - reference `kraken.exe`/`krakenc.exe` binaries
- `Printf`, `LinearAlgebra` - used by the `test/reference/` harness

> **Stdlibs must be listed too, and `Pkg.test()` is the only check that catches a missing one.**
> `--project=test` resolves against a manifest that already contains every stdlib as an indirect
> dependency, so `using LinearAlgebra` there succeeds even when `test/Project.toml` never declared
> it. `Pkg.test()` builds the environment from the declared dependencies only and fails with
> `ArgumentError: Package LinearAlgebra not found in current path`. Run the full `Pkg.test()` before
> pushing; a green single-file run proves less than it appears to.
