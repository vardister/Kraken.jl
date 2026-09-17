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

There is no single M4/M5 figure in that table on purpose: since Milestone 4 the whole suite in one
call exceeds every limit the MCP imposes (see the bottom of this file), so it is measured per file
instead. CI still runs it as one `Pkg.test()` — that path has no such limit.

Per file, so a silent drop in coverage is visible in a diff. The M5 column is after task 5.3, with an
Acoustics Toolbox checkout present, and was measured **per file** — the whole suite in one call is not
obtainable through the MCP (see "Running the suite through the kaimon MCP" below):

| File | End of M1 | End of M2 | End of M3 | End of M5.5 |
|---|---|---|---|---|
| `environment_tests.jl` | 39 | 161 | 161 | 215 |
| `integration_tests.jl` | 98 | 98 | 98 | 134 |
| `numerical_methods_tests.jl` | 73 | 73 | 73 | 96 |
| `automatic_differentiation_tests.jl` | 48 | 48 | 48 | 48 |
| `reverse_ad_tests.jl` | — | — | — | 359 |
| `fortran_reference_tests.jl` | — | — | 545 | 793 |
| `performance_tests.jl` (opt-in) | 24 | 24 | 24 | 24 |

**1645 in total at the end of Milestone 5.5**, all green, measured file by file with an Acoustics
Toolbox checkout present. (At the end of 5.3 it was 1558; 5.4 added 50 attenuation-AD assertions and
5.5 added 37 for the lossy standard environments.)

**Counting the TestItems files from a worktree needs a path filter.** `@run_package_tests` walks the
whole package directory, and `.claude/worktrees/` is inside it — so a filter that only matches on
`endswith(t.filename, "integration_tests.jl")` silently runs *every* sibling worktree's copy too and
reports their sum. Add the worktree name:

```julia
@run_package_tests filter = t -> occursin("plan-5-attenuation", t.filename) &&
                                 endswith(t.filename, "integration_tests.jl")
```

Without it this table read 232 rather than 134 for `integration_tests.jl`, which looks like coverage
that is not there.

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

Coverage as of 2026-09-16, after Milestone 6.4 — **211 of 402 `.env` files parse, up from 167** after
Milestone 5.1. Milestone 6 lifted the rigid and vacuum boundaries and the n²-linear and cubic-spline
interpolators:

| Blocker | after 5.1 | after 6.4 | Status |
|---|---|---|---|
| **parses** | **167** | **211** | |
| top boundary | 65 | 65 | all `A` (acousto-elastic half-space above) — out of scope, plan task 6.1 |
| bottom half-space is not the fastest medium | 31 | 34 | leaky modes (M5 stretch) |
| added volume attenuation (THORP / Francois-Garrison / biological) | 27 | 27 | `TopOpt(4:4)`, a separate feature |
| profile does not start at the surface | 1 | 21 | — |
| elastic layer | 7 | 9 | out of scope |
| interfacial roughness | 4 | 6 | out of scope |
| bottom boundary | 50 | 6 | all `F` (reflection-coefficient file), which `kraken.f90` itself rejects |
| SSP interpolation over a varying profile | 28 | 1 | a cubic spline over several media (`Munk/MunkB_ray.env`) |
| power-law attenuation (`TopOpt(3:3) == 'm'`) | 1 | 1 | needs per-medium β and f_T records |
| analytic SSP | 1 | 1 | — |
| not a KRAKEN deck (BELLHOP3D `'H'`/`'Q'` SSP options, malformed) | 20 | 20 | n/a |

Some blockers *rose*. That is files the reader used to reject early, for their bottom option or
interpolator, now getting far enough to be judged on something else — mostly BELLHOP decks whose
profile starts below the surface.

**Parsing is not reproducing**, so the milestone was also measured the way its success criterion is
worded. Of the files that parse, the ones `kraken.exe` itself can run *as shipped* were solved by both
and counted as reproduced when every compared `kᵣ` agrees within 1e-3 (the tolerance the lossy cases
use). Most unrunnable decks are BELLHOP inputs whose NMESH is "too coarse" for KRAKEN:

| | parses | parses and `kraken.exe` runs it | Kraken.jl reproduces `kᵣ` |
|---|---|---|---|
| after 5.1 (reader on `revival`, 2026-09-16) | 167 | 81 | **81** |
| after 6.4 | 211 | 100 | **97** |

The three new decks that do not reproduce are all perfect bottoms, and are the two defects described
under "Boundary conditions and SSP interpolation validated against Fortran" below: `wedge/wedge.env`
(vacuum, 4.2e-3 on its grazing modes) and `Dickins/Precalc/DickinsK.env` and
`Bellhop3DTests/DoubleSeamount/DoubleSeamount3D_ray.env` (the mode-cutoff crash). Mode-shape
correlation was deliberately left out of that count: at a thousand modes the writer's receiver grid,
capped at 2001 points, cannot resolve the mode shapes, so `MunkS_500Hz.env` and the `DoubleSeamount`
decks read ~0 on both branches.

Regenerate this with `KrakenReference.categorize_env_tree` and `print_env_tree_report`. A file that
uses an unsupported feature is *named*, never approximated — the whole point is that a case Kraken.jl
cannot model fails with "unsupported feature: top boundary (acousto-elastic halfspace)" rather than
silently mis-parsing into a plausible environment that then disagrees with Fortran for reasons nobody
can find.

Three caveats the suite encodes rather than hides:

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
- **A record ends where Fortran stops reading it, not at end of line.** `READ( ENVFile, * ) BotOpt,
  Sigma` is a list-directed read of exactly two items and **discards the rest of the record**. The
  reader used to take the *last* number on the line as SIGMA, so `SedAtten/calibS_0.6dB.env`'s
  `'A'  0.0 2.5 2000` — whose trailing pair only a `'m'` attenuation unit would go back for — read as
  2 km of interfacial roughness and the file was rejected. Fixed in 5.1 and pinned by
  `"M5.1: the bottom-option record ends after SIGMA"`. The same rule applies to any record where the
  toolbox files carry optional trailing parameters.

### Attenuation validated against Fortran (plan task 5.3)

Milestone 5 makes `kraken_jl` return **complex** wavenumbers when an environment declares
attenuation. `Im(kᵣ)` is compared separately from `Re(kᵣ)` and gets its own, looser tolerance —
`KrakenReference.max_alpha_reldiff` alongside `max_kr_reldiff`.

Max relative difference over all modes, measured 2026-08-09. The first four rows are `pekeris_env` at
100 Hz with attenuation added; the rest are read from the toolbox tree.

| case | modes | max rel Δ Re(kᵣ) | max rel Δ Im(kᵣ) | tolerance asserted |
|---|---|---|---|---|
| `pekeris_env`, lossless (control) | 5 | 1.7e-9 | 0 exactly, both sides | 1e-6 / 1e-12 |
| `pekeris_env` + 0.5 dB/λ in the water | 5 | 1.3e-4 | 4.0e-3 | 1e-3 / 2e-2 |
| `pekeris_env` + 0.05 dB/λ in the half-space | 5 | 1.3e-5 | 7.9e-4 | 1e-4 / 5e-3 |
| `pekeris_env` + 0.5 dB/λ in the half-space | 5 | 5.5e-4 | 1.0e-1 | 3e-3 / 3e-1 |
| `SedAtten/calibS_noloss.env` (control) | 11 | 6.0e-10 | 0 exactly, both sides | 1e-6 / 1e-12 |
| `TLslices/atten.env`, 10 Hz, 0.001 dB/(km·Hz) | 44 | 7.2e-7 | 1.8e-3 | 1e-5 / 1e-2 |
| `SedAtten/calibK.env`, 250 Hz, 0.5 dB/λ | 11 | 1.8e-4 | 1.1e-2 | 1e-3 / 5e-2 |
| `SedAtten/calibS_0.6dB.env`, 250 Hz, 0.6 dB/λ | 11 | 2.3e-4 | 1.5e-2 | 1e-3 / 8e-2 |

Three things about that spread are worth understanding before changing any of these numbers.

- **`Im(kᵣ)` is a harder quantity than `Re(kᵣ)` on both sides.** It is the small part of the complex
  number; the only usable Fortran source for it is the single-precision `.mod` (the `.prt` prints it
  with `G10.2`, i.e. two digits); and *neither* solver Richardson-extrapolates it — both evaluate the
  perturbation on their own coarsest mesh, which is what `kraken.f90` does by calling `Vector` only
  when `iSet == 1`.
- **The two solvers agree to first order in α exactly.** Shrinking the half-space attenuation drives
  the disagreement down as α², not α: 7.1e-3 at 0.5 dB/λ → 6.9e-5 at 0.05 → 1.9e-6 at 0.005, where it
  hits the `.mod`'s single-precision floor. A formula error would leave a floor proportional to α.
  `"M5.3: the two solvers agree to first order in α"` asserts exactly this, and it is the strongest
  correctness statement in this section — much stronger than any single tolerance.
- **What is left is second order, and it is a limit of the method, not of the implementation.** The
  half-space term goes as `Im √(γ² + 2iω α_b/c_b)`, whose expansion parameter is
  `v = 2ω α_b/(c_b γ²)` with `γ² = kᵣ² − (ω/c_b)²`. `γ → 0` at the bottom cutoff, so `v` is large
  exactly where the mode is least trapped: it is 1.5e-1 for `calibK` mode 1 and **3.2** for
  `pekeris` mode 5. Both codes are first-order methods keeping different second-order terms, so at
  `v ≈ 3` neither is more right — the full complex solve (`krakenc.exe`, plan task 5.6) is the tool
  for that regime. The volume term has no `γ` in it, which is why moving the same 0.5 dB/λ from the
  half-space into the water improves agreement by a factor of 12 on the same waveguide.

#### A second, unrelated limit: quadrature at a discontinuity (task 5.5)

The lossy *standard environments* added in 5.5 are in the suite too. One of them — the only case with
a lossy *interior* medium — exposed an accuracy limit unrelated to the half-space cutoff above, which
task 5.7 then fixed. Both columns are kept because the contrast is the evidence:

| environment | max rel Δ Re(kᵣ) | Im(kᵣ) before 5.7 | Im(kᵣ) after 5.7 | tolerance asserted |
|---|---|---|---|---|
| `pekeris_env(; α0=0.2)` | 2.1e-5 | 7.9e-4 | 8.1e-4 | 1e-4 / 5e-3 |
| `pekeris_env(; αb=0.2)` | 1.6e-4 | 5.8e-2 | 5.8e-2 | 1e-3 / 2e-1 |
| `one_layer_env(; α1=0.4)` | 3.0e-5 | **8.1e-2** | **2.9e-3** | 1e-4 / 1e-2 |
| `one_layer_env(; α0=0.1, α1=0.4, αb=0.1)` | 4.3e-5 | 8.7e-3 | 7.5e-4 | 1e-4 / 5e-3 |

Exactly the pattern the fix predicts: the two `pekeris` rows are a single water medium over a
half-space, have no interior interface to straddle, and are untouched — their residual is the
bottom-cutoff limit above. The two `one_layer` rows improve 27x and 12x.

**What was wrong.** `ρ` jumps at every interface, and `α` jumps wherever the loss is confined to a
layer, so the integrands of *both* the energy normalization and the attenuation perturbation are
discontinuous. Both ran a single `integral_trapz` over the whole flattened mesh, and one straddled
interval per interface is enough to drop the entire quadrature from second order to first. The
resulting error pattern was **inverted** from the half-space cases — worst on the *best*-trapped
mode, whose loss comes entirely from an exponentially small tail inside the sediment, sampled exactly
where the quadrature was weakest.

**The fix.** `normalize_mode` and `modal_attenuation` integrate medium by medium, as `kraken.f90`'s
`Normalize` does. `medium_mesh` / `medium_mode` / `medium_property` in `src/kraken_core.jl` extend
each medium's samples up to its own top interface: the mode is continuous there so its value carries
over from the medium above, while the *profiles* are extrapolated from inside the medium as
`2p₁ - p₂` — exact, because the profiles are piecewise linear and the augmented mesh is uniform. The
two functions had to move together, since the perturbation is a ratio against the normalization and
they must stay weighted alike.

**Observed convergence order on `one_layer_env(; α1=0.4)`, mode 1:**

| | 20→40 | 40→80 | 80→160 | 160→320 |
|---|---|---|---|---|
| before 5.7 | 1.02 | 1.07 | 1.11 | 1.29 |
| after 5.7 | 2.00 | 2.00 | 2.02 | 2.07 |

**Measure that order against Kraken.jl's own fine solution, not against `kraken.exe`.** On this case
`kraken.exe`'s automatic mesh carries about **1.6e-3** of its own discretization error, and that is a
floor rather than a slope: measured against Fortran the apparent order collapses to 0.57, 0.19, 0.05
as our error drops below its. `"M5.7: the perturbation integral converges at second order"` uses a
self-reference for exactly this reason, and says so in a comment.

`VolAtt`, named in the plan, is deliberately **not** in the table: every file in it puts an
acousto-elastic half-space *above* the surface and gives the bottom the water's own sound speed, so
there is no trapped spectrum to compare — they are free-space TL cases, not modal ones, and two of
them additionally use `TopOpt(4:4)` volume-attenuation laws. `TLslices/atten.env` takes their place
and is a better test anyway: 44 modes, loss in *both* media, and the only case exercising dB/(km·Hz).

### Boundary conditions and SSP interpolation validated against Fortran (plan task 6.5)

Each option Milestone 6 added is compared against `kraken.exe` in `fortran_reference_tests.jl`.
Measured 2026-09-16.

**Boundary conditions** — `pekeris_env`, which isolates the boundary rows because the column is
isovelocity:

| top / bottom | Hz | modes | max rel Δkᵣ | min corr | tolerance asserted |
|---|---|---|---|---|---|
| rigid / half-space | 100 | 5 | 6.2e-10 | 0.9999999 | 1e-8 / 0.9999 |
| vacuum / rigid | 100 | 13 | 1.4e-10 | 0.9999944 | 1e-8 / 0.9999 |
| vacuum / vacuum | 100 | 13 | 6.3e-9 | 0.9999935 | 1e-7 / 0.9999 |
| rigid / rigid | 100 | 14 | 3.3e-10 | 0.9999928 | 1e-8 / 0.9999 |
| rigid / vacuum | 100 | 13 | 1.4e-10 | 0.9999939 | 1e-8 / 0.9999 |
| vacuum / vacuum | 50 | 6 | 2.2e-10 | 0.9999932 | 1e-8 / 0.9999 |
| vacuum / rigid | 50 | 7 | 6.1e-5 | 0.9999905 | 2e-4 / 0.9999 |
| `munk_env`, n²-linear, vacuum / rigid | 10 | 66 | 7.9e-6 | 0.9999854 | 5e-5 / 0.9999 |

The 50 Hz rigid-bottom outlier is mode 7, which grazes at `kᵣ = 0.046` against `ω/c = 0.209`. A
small absolute error in `kᵣ²` reads large relative to a small `kᵣ`.

**SSP interpolation** — `munk_env` cannot test this. Its 100 m sampling of a smooth profile barely
separates the interpolators: at 25 Hz the n²-linear solve matched Fortran's n²-linear run to 6.4e-6 and
its *C-linear* run to 6.8e-6. The suite uses a coarse duct instead, five samples 50 m apart
(1540, 1500, 1480, 1500, 1530 m/s over a 1700 m/s half-space), and checks the whole matrix. Max rel
Δkᵣ at 50 Hz, Julia interpolation (rows) against Fortran interpolation (columns):

| | Fortran `C` | Fortran `N` | Fortran `S` |
|---|---|---|---|
| `:c_linear` | **7.7e-8** | 1.2e-4 | 1.8e-3 |
| `:n2_linear` | 1.2e-4 | **7.0e-8** | 1.8e-3 |
| `:cubic_spline` | 1.8e-3 | 1.8e-3 | **5.2e-9** |

C-linear and n²-linear are told apart by three orders of magnitude, so a wrong interpolator fails.
100 Hz is within 4× of every entry (the spline diagonal is 9.6e-9 there). Asserted: diagonal < 1e-6 for
all three, off-diagonal > 1e-5 for C and N and > 1e-3 for the spline. The spline diagonal was 3.1e-4
until task 6.8 — see below.

**The toolbox's newly readable decks** (CLOW/CHIGH from the file; asserted in
`"M6.5: the toolbox's newly readable options against kraken.exe"`):

| deck | option | Hz | modes (Fortran) | max rel Δkᵣ | min corr | tolerance asserted |
|---|---|---|---|---|---|---|
| `TLslices/pekeris.env` | N | 10 | 44 | 1.8e-7 | 0.9999930 | 1e-6 / 0.9999 |
| `Noise/Pekeris/pekeris.env` | N | 300 | 26 | 3.7e-5 | 0.9999939 | 1e-4 / 0.9999 |
| `MunkLeaky/MunkK1525.env` | N | 50 | 28 | 3.2e-6 | 0.9999992 | 1e-5 / 0.9999 |
| `Gulf/gulf_rd.env` | N | 50 | 63 | 1.1e-6 | 0.9999929 | 1e-5 / 0.9999 |
| `Munk/MunkK.env` | N | 50 | 102 | 6.1e-5 | 0.9943636 | 2e-4 / 0.99 |
| `Munk/MunkS.env` | S | 50 | 102 | 6.1e-5 | 0.9944935 | 2e-4 / 0.99 |
| `wedge/wedge.env` | vacuum bottom | 25 | 59 | 4.2e-3 | 0.9998710 | 1e-2 / 0.999 |

The Munk decks' worst modes are 99–102, the last trapped modes at the half-space cutoff, and **it is
not the interpolator**: the same decks forced to C-linear on both sides give the same 6.1e-5 and
0.9944. `wedge.env`'s worst modes are its grazing ones near `kᵣ = 0` (mode 59: `kᵣ = 0.0147` at
`ω/c = 0.105`). `MunkK1525` and `gulf_rd` narrow CHIGH to 1525 m/s, so Fortran reports only their
slow modes while Kraken.jl finds every trapped one; the leading modes are what is compared.

#### Two defects this found

Both were in `src/`, outside a test task's reach, and became plan tasks 6.7 and 6.8. Both are fixed,
and the `@test_broken` that pinned them are now plain `@test`.

- **Fixed in 6.7 — a perfect bottom could crash the solve when a mode sat at `kᵣ = 0`.** `richard_extrap` takes the
  square root of an extrapolated `kᵣ²` that has crossed zero (`DomainError`); KRAKEN discards such a
  mode (`kraken.f90` keeps only `kᵣ² > ω²/cHigh²`). For the 100 m Pekeris column the cutoffs are the
  multiples of `c/2D = 7.5 Hz` with a vacuum bottom and the odd multiples of 3.75 Hz with a rigid one.
  Swept over 20–200 Hz in 0.25 Hz steps: 10 of 721 frequencies fail for vacuum (75, 90, 120, 127.5,
  … Hz — 10 of the 24 cutoffs in range) and 6 for rigid (93.75, 168.75 and 198.75 Hz, plus a
  `SingularException` in inverse iteration's LU at 183.25–183.75 Hz, next to the 183.75 Hz cutoff).
  Every failure is at or beside a cutoff, and round frequencies are the likely ones to be asked for. `Dickins/Precalc/DickinsK.env` (rigid, 230 Hz)
  and `Bellhop3DTests/DoubleSeamount/DoubleSeamount3D_ray.env` hit it in the toolbox tree.

  After 6.7 the same sweep has **no failures** for either bottom, and at all 14 former failure
  frequencies Kraken.jl and `kraken.exe` agree on the mode count exactly and on `kᵣ` to ≤ 3.6e-9. A
  mode whose extrapolated `kᵣ²` is below `√eps · max(ω/c)²` is dropped, as KRAKEN drops
  `kᵣ² ≤ ω²/cHigh²`. The singular LU was a shift sized to `kᵣ` rather than to the scaled diagonal; the
  retry that replaces it runs only when the LU fails, and 44 pre-change solutions (every standard
  environment lossless and lossy, every boundary pair, both new interpolators) are bit-identical.
  `DickinsK.env` now solves (1226 Julia modes; `kraken.exe` reports 1212 because the deck's
  `CHIGH = 10000` cuts off `kᵣ < 0.1445`). Its modes 1–264, trapped above the sediment, agree to
  3e-8; beyond, they reach into 1000 m of 0.5 dB/λ sediment and drift to 1.6e-3 in `kᵣ`. That is the
  lossy-medium limit of 5.3/5.5, not this: with the attenuation zeroed, all 1212 agree to 2.6e-6.
- **Fixed in 6.8 — `:cubic_spline` was a different spline from KRAKEN's.** Kraken.jl used DataInterpolations'
  *natural* spline; `cCubic` calls `CSPLINE` with `IBCBEG = IBCEND = 0`, which `misc/splinec.f90`
  documents as *not-a-knot*. On the duct that is the whole 3.1e-4. Sampling each end condition's
  spline at 0.5 m and solving it as a C-linear profile, the not-a-knot one matches Fortran's spline
  run to 2.1e-7 and the natural one matches Kraken.jl's to 2.0e-7. `"the spline gap is the end
  condition"` asserted both. On a finely sampled profile the two end conditions barely differ
  (`Munk/MunkS.env`, splined, agreed to 5.4e-5; its C-linear control to 6.1e-5).

  After 6.8 the spline is `NotAKnotSpline`, a line-for-line transcription of `CSPLINE` with
  `IBCBEG = IBCEND = 0` (two points give the line, three the parabola). The duct's spline diagonal is
  5.2e-9 at 50 Hz and 9.6e-9 at 100 Hz, and the sampled not-a-knot profile now matches Kraken.jl's own
  spline solve as well as Fortran's. `MunkS.env` reads 6.1e-5 / 0.99449, the same as its C-linear
  control, so what remains there is the cutoff modes, not the interpolator.

### AD through a lossy solve (plan task 5.4)

`modal_attenuation` needed no rule of its own — it is traced arithmetic downstream of the two
functions that carry rules — and Zygote agrees with ForwardDiff to ~1e-11 on the full parameter
gradient of a lossy Pekeris solve, `∂/∂α` included. `"Reverse-mode AD through attenuation"` in
`test/reverse_ad_tests.jl` is the 50 assertions that make that checkable. Two limits it pins:

- **At `α = 0` the two AD modes return opposite one-sided derivatives, silently.** `imag(kᵣ)` is
  identically zero below zero attenuation and linear above, so the point is a kink. ForwardDiff's
  `iszero` on a `Dual` sees the seed, `is_lossy` takes the lossy branch, and it reports the
  right-hand derivative (-0.0372 on the Pekeris case). Zygote evaluates `is_lossy` on the primal,
  gets `false`, the attenuation path never reaches the tape, and it reports the left-hand one, zero.
  One hair above zero they agree to 1e-15. **Differentiate at a nonzero attenuation.** If you are
  chasing a "broken rule" because two backends disagree on an attenuation derivative, check this
  first.
- **Mooncake cannot trace the complex path.** `add_attenuation`'s `sqrt` of a `Complex` trips
  `ArgumentError: It is not permissible to bitcast to a differentiable type during AD`. The failing
  call is inside Mooncake's own `Complex` handling, so no rule in `src/kraken_ad.jl` fixes it —
  reverse mode over a lossy environment means Zygote today. Pinned as a *specific* expected failure
  matching on `"bitcast"`, so Mooncake gaining complex support shows up as that test going green.

The convention, for anything downstream: **real parameters in, real loss out, complex only in
between.** `sum(imag, kr)`, `sum(real, kr)` and `sum(abs2, kr)` all work on every backend. A
complex-valued loss is refused by Zygote rather than resolved by picking a conjugation convention,
and that refusal is asserted.

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

### Running the suite through the kaimon MCP

`CLAUDE.md` directs every Julia invocation through the `kaimon` MCP rather than Bash, so that the run
is visible in the shared REPL. That works for everything here *except* the suite as a single call.
Since Milestone 4 the full run is past six minutes, and both whole-suite routes fail at that length
(measured 2026-08-09):

| Route | What happens |
|---|---|
| `run_tests(project_path=…)` | The MCP **transport drops mid-call** — "response for tool run_tests was lost". Three times out of three. The server itself stays up; `ping` answers normally straight afterwards. |
| `ex(e="using Pkg; Pkg.test()")` | Killed by the gate's **10-minute inactivity timeout**, at exactly 10m00s. `ex` strips stdout, so the gate sees no output and concludes the session is stuck. |

Neither is the suite hanging. Run the files individually instead, from a session rooted at `test/` —
each `include` returns its own `Test Summary`, and none of them is long enough to trip either limit:

```julia
include(joinpath(@__DIR__, "fortran_reference_tests.jl"))          # ~1-3 min
include(joinpath(@__DIR__, "numerical_methods_tests.jl"))          # ~3 s
include(joinpath(@__DIR__, "automatic_differentiation_tests.jl"))  # ~2.5 min
include(joinpath(@__DIR__, "reverse_ad_tests.jl"))                 # ~5 min
using TestItemRunner
@run_package_tests filter = t -> endswith(t.filename, "environment_tests.jl")
```

A **subprocess launched from inside `ex`**, with its output read back as a value, is how to get a
`Test Summary` out of a run at all — `ex` strips stdout, so a `Test Summary` printed by an `include`
goes into the void, and a TestItems file run with `@run_package_tests` otherwise reports nothing:

```julia
cmd = `$(Base.julia_cmd()) --project=test -e "using TestItemRunner; @run_package_tests filter=t->endswith(t.filename, \"integration_tests.jl\")"`
out = read(pipeline(ignorestatus(cmd), stderr=devnull), String)
filter(l -> occursin("Test Summary", l) || occursin("Fail", l), split(out, '\n'))
```

**It does not lift the 10-minute limit** — tried on the whole suite on 2026-08-09 and killed at
exactly 10m00s, same as a direct `Pkg.test()`. Capturing the subprocess's stdout with `read` is
precisely what keeps the gate from seeing any activity. So the per-file rule above stands: use the
subprocess to *see* results, and still run one file at a time. `ex` promotes anything past 30 s to a
background job; poll `check_eval` sparingly.

Two things to watch. A TestItems filter that matches nothing still reports green, so confirm it
selected something before believing the result — `filter = t -> (push!(items, t.filename); false)`
counts the items it would have run (22 across `environment_tests.jl` and `integration_tests.jl`).
And per-file runs are weaker evidence than `Pkg.test()`: they resolve against `test/Manifest.toml`,
so they cannot catch an undeclared dependency (see the stdlib note below). Say which one you ran.

Wall-clock times above vary by a factor of three depending on how many Julia sessions are live —
the same `fortran_reference_tests.jl` took 1m07 on an idle machine and 3m03 with three sessions and
background jobs competing. Shut down sessions you are not using (`manage_repl(command="shutdown")`)
before reading anything into a timing.

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
