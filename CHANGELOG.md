# Changelog

All notable changes to Kraken.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses
[semantic versioning](https://semver.org/) — in 0.x, a minor bump may be breaking.

## [Unreleased]

### Changed

- **`kraken_jl` is 4.9 to 11.8 times faster than in 0.4.0, and faster than Fortran `kraken.exe` on
  every benchmark case.** Munk at 50 Hz (102 modes) went from 355 ms to 30.1 ms, against 31.4 ms for
  a local gfortran build of `kraken.exe` with its process launch subtracted; Pekeris at 500 Hz from
  12.9 ms to 1.31 ms, against 3.6 ms. Allocation per Munk solve fell from 67.9 to 10.2 MiB.
  `test/benchmark_vs_fortran.jl` measures all of this. The changes:
  - **Inverse iteration** tests convergence up to sign. A shift on the far side of the eigenvalue
    flips an already-converged vector on every solve, and the old componentwise test never passed:
    33 of Munk's 102 modes ran the full 50 iterations. The matrix is factored once per mode, by a
    `gttrf`-style LU with reciprocal pivots into reused task-local buffers, and the convergence test
    takes one sweep instead of three.
  - **The default root solver is `SafeguardedNewton()`**: Newton's method on the Sturm determinant
    inside the bisection bracket, bisecting whenever a step would leave it. On the first Munk mesh
    it takes 7.6 determinant evaluations a mode, against 10.4 for `Brent()` and about 30 for `ITP()`,
    whose steps on these brackets are almost all bisections. Each finer mesh starts from the
    previous mesh's roots (`find_kr`'s new `guesses` keyword), and then takes 3.9.
  - **Normalization** is a dot product with weights computed once per mesh (`mode_energy_weights`),
    instead of resampling the density and rebuilding the medium meshes for every mode: 107 to 9.6
    µs a mode.
  - **`det_sturm` and `sturm_and_slope`** multiply by the rescaling factor only when it is not 1, and
    handle the bottom row outside the loop. Both results are bit for bit the same; on a 6250-point
    mesh `det_sturm` went from 41 to 18.7 µs and `sturm_and_slope` from 38.3 to 21.9 µs.
  - `kraken_jl` no longer copies the mode matrix twice on return.
- **`kraken_jl`'s default `abstol` and `reltol` are `1e-10`, not `1e-6`.** They apply only to a
  NonlinearSolve `method`; `SafeguardedNewton()` converges to 16 ulps. With `Brent()` at `1e-6`,
  Munk's wavenumbers were 1.3e-6 from `kraken.exe`; at `1e-10`, 3.7e-8. `compare_with_fortran`'s
  defaults follow.
- **The solver's branches compare plain values**, through `primal_value`, so a ForwardDiff solve
  makes the same decisions as the plain one. ForwardDiff's `<` breaks ties between equal values by
  the partials; on Pekeris at 184 Hz that sent the two down different paths and mode 9 came out
  3e-15 apart.
- **`SafeguardedNewton()` also stops once Newton has stalled at the determinant's rounding noise.**
  On SWellEx-96's 32032-point finest mesh at 500 Hz the steps settled between 28 and 50 ulps, never
  met the 16-ulp test, and every mode ran 200 iterations: the solve took 9.34 s, now 0.75 s.

### Added

- `test/solver_internals_tests.jl`: unit tests for the pieces above, each against an independent
  reference (270 assertions, about 10 s).
- The SWellEx-96 environment at 500 Hz (432 modes) in `test/benchmark_vs_fortran.jl`, from the
  experiment's published geoacoustic model and CTD cast `i9606`.

### Known issues (found with SWellEx-96, 2026-10-16)

- **Mesh sizing.** Every layer's mesh spacing is `c_b / (20 f)`, from the half-space speed. With a
  5200 m/s half-space under 1488 m/s water that is about 6 points per water wavelength; `kraken.exe`
  sizes each medium from its own sound speed and refuses this mesh as too coarse. Wavenumbers still
  agree to 7.7e-5, but Kraken.jl needs 5 refinement levels where `kraken.exe` needs 3, and is 2.9
  times slower on this case (0.75 s against 0.26 s).
- **Attenuation near the half-space cutoff.** `Im(kᵣ)` differs from `kraken.exe` by a median of 7%
  on the default mesh and 4% at 20 points per water wavelength, and by up to about 80% for modes
  with phase speeds near the 5200 m/s cutoff. The cause is not yet identified.

## [0.4.0] — 2026-10-16

Everything since 0.3.0 is new capability: reverse-mode AD, attenuation, the boundary conditions and
profile interpolations KRAKEN offers, the pressure field, a documentation site, and plotting. Each is
cross-validated against unmodified Fortran KRAKEN (`kraken.exe`, `field.exe`) in the test suite.

### Breaking

- **A lossy environment now returns complex wavenumbers.** 0.3.0 read the attenuation columns of
  `ssp` and `sspHS` and ignored them. They are now used, and `kraken_jl(env, f).kr` is `ComplexF64`
  whenever any attenuation is nonzero, with `Im(kᵣ) < 0`. A lossless environment still returns real
  `Float64` wavenumbers, bit-identical to 0.3.0.
- **An interface written as two SSP rows at the same depth now solves correctly**, so its results
  change. The water column's last mesh point used to take the sediment's sound speed and density,
  which put `kᵣ` up to 8% from `kraken.exe` on a shallow shelf. The `z0`, `z0 + eps(z0)` spelling
  that the standard environments use was always right and is unchanged.
- `dev/kraken_broadband.jl` is gone. Broadband synthesis is package API now:
  `broadband_field` and `synthesize_pulse`.

### Added

- **Reverse-mode AD** through Zygote.jl and Mooncake.jl (the latter via a package extension). The
  root finder and the eigenvectors carry implicit-function rules, and the sound-speed and density
  interpolants carry rules that include their knot depths. A gradient with respect to a whole
  sound-speed profile costs a fixed multiple of one solve, however many points the profile has.
- **Attenuation**, added by first-order perturbation as `kraken.exe` does, in all six of KRAKEN's
  unit conventions (`atten_units`). The standard environments take `α` keywords.
- **Boundary conditions:** `PressureRelease`, `RigidBoundary` and `AcousticHalfspace`, through
  `UnderwaterEnv`'s `top_bc` and `bottom_bc`.
- **Sound-speed interpolation:** `:c_linear` (default), `:n2_linear` and a not-a-knot
  `:cubic_spline` matching KRAKEN's, through `ssp_interp`.
- **The pressure field:** `acoustic_field`, `transmission_loss` (coherent and incoherent) and
  `mode_amplitudes`, in `field.exe`'s conventions, plus broadband synthesis. All of it is
  differentiable.
- **Plotting**, through the Makie extension: `plot_tl`, `plot_tl_range`, `plot_dispersion`,
  `plot_arrivals` and `plot_ssp!`, alongside `plot_modes` and `plot_ssp`.
- **A documentation site**, with a getting-started guide, the environment format, AD, a validation
  page whose numbers the docs build measures by running Fortran KRAKEN, and an API reference in which
  every exported name has a docstring, and every example outside the plotting functions is a doctest.
- `examples/`, with its own environment: deep-water modes, dispersion on a lossy shelf, a threaded
  frequency sweep, and a tour of reverse-mode AD. CI runs them all, and runs every code block in the
  README.

### Fixed

- **Forward-mode derivatives of the wavenumbers were slightly wrong.** ForwardDiff went through the
  root finder's iterations, so a wavenumber's derivative came from the last bracketing step, and
  group speeds were off by up to ~0.1% at scattered frequencies. A new ForwardDiff extension
  (ForwardDiff is a weak dependency) applies the same implicit-function derivative as reverse
  mode, and the two modes now agree to rounding.
- Reverse mode differentiated the wrong side of a same-depth interface, which put some gradients
  15% off. It is fixed by the interface change above.

### Known limitations

- Shear waves: `cs` and `αs` are read and ignored.
- Leaky modes: attenuation is perturbational. The complex solve of `krakenc.exe` is not implemented.
- Range dependence.
- `:cubic_spline` over more than one medium. KRAKEN splines each medium separately, while Kraken.jl
  runs one spline through the whole column.

## [0.3.0] — 2026-08-28

The first release of the revived package. 0.2.3 shipped code whose test suite could not run, whose
docs build had never executed, and which pulled in a dozen dependencies it did not use. This release
is that tree repaired, pruned, put under continuous integration, and cross-validated against
unmodified Fortran KRAKEN.

### Breaking

- The Fortran-backed `kraken` function and the `EnvKRAKEN` API are gone. The entry point is
  `kraken_jl(env, freq)`, which returns a `NormalModeSolution` (`kr`, `modes`, `env`, `props`).
  Nothing in the package links against Fortran code or ships a shared library.
- Plotting moved out of the package into the `KrakenMakieExt` package extension. `plot_modes` and
  `plot_ssp` now require `using CairoMakie` (or `GLMakie`); without a backend loaded they raise an
  actionable error rather than a `MethodError`. `src/kraken_plots.jl` was deleted — it consumed the
  `Dict` returned by the removed Fortran path and could not have run against the current API.
- The standard-environment helpers were renamed: `pekeris_env`, `one_layer_env`, `munk_env`, and
  friends.
- Ten unused dependencies were dropped — CairoMakie, DSP, Distributed, FFTW, FiniteDiff, JLD2,
  NamedArrays, Parameters, ProgressMeter, and QuadGK — and Makie became a weak dependency. Load and
  install time drop accordingly.
- Broadband pulse synthesis moved to `dev/kraken_broadband.jl`. It is a script you `include` by
  hand, not package API.

### Added

- **Cross-validation against unmodified Fortran KRAKEN**, on every push. The binaries come from
  [`AcousticsToolbox_jll`](https://github.com/JuliaBinaryWrappers/AcousticsToolbox_jll.jl), so no
  Fortran toolchain is needed — including in CI — and the comparison is driven over the toolbox's
  own `.env`/`.mod`/`.prt` file interface. Across the five standard environments at 25–400 Hz, the
  largest relative wavenumber difference is **2.7e-5** and the smallest mode-shape correlation is
  **0.99995**. Of the 402 environment files shipped with the Acoustics Toolbox, the 65 that use only
  features Kraken.jl models today all agree.
- Acoustics Toolbox `.env` files are both **read and written** (`test/reference/`), so existing
  environment files can be reused directly.
- GitHub Actions CI: the test suite on Julia 1.10 / lts / 1 across Ubuntu and macOS, a
  JuliaFormatter gate, and a documentation build.
- An opt-in performance suite behind `KRAKEN_RUN_PERFORMANCE_TESTS=true`.

### Fixed

- The test environment was wired such that `Pkg.test()` could not run at all. It runs, and it is
  green.
- Five real bugs, each pinned by a regression test verified to fail on the pre-fix source:
  - `Base.show(::SampledDensity1D)` referenced a nonexistent `.type` field, so displaying any
    density profile — or any environment — threw a `FieldError`.
  - `pressure_f` called a five-argument `get_modal_function` that does not exist; every call raised
    a `MethodError`.
  - `det_sturm`'s `stop_at_k` early exit was an expression statement rather than a `return`, so the
    option ran the loop to completion and silently did nothing. It was removed, along with a
    documented `return_det` keyword that never existed.
  - The two `UnderwaterEnv` constructors read the environment depth from different places
    (`layers[end, 3]` versus `ssp[end, 1]`) and disagreed whenever the SSP was sampled past the last
    layer boundary. Both now use `layers[end, 3]`.
  - `munk_env` returned its SSP as a lazy `Transpose`, which made `UnderwaterEnvFORTRAN(munk_env()...)`
    error outright.
- The documentation build, broken for the entire life of the repository (`docs/make.jl` said
  `using KRAKEN`; the module is `Kraken`) and never executed by anything. It now runs in CI.
- Two incorrect assertions in the Munk-profile tests.
- Tracked binaries and generated run artifacts (`.mod`, `.prt`, `.shd`, `.dylib`, …) were removed
  from the repository and ignored.

### Known limitations

Not yet modelled — see the README for the current list:

- Compressional wave attenuation (the single largest gap against the Acoustics Toolbox test cases)
- Shear wave properties
- Boundary conditions other than a pressure-release surface over a fluid half-space
- SSP interpolation other than C-linear
- Reverse-mode automatic differentiation (forward mode, via ForwardDiff.jl, works today)

[0.4.0]: https://github.com/vardister/Kraken.jl/releases/tag/v0.4.0
[0.3.0]: https://github.com/vardister/Kraken.jl/releases/tag/v0.3.0
