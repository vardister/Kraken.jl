# Changelog

All notable changes to Kraken.jl are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project uses
[semantic versioning](https://semver.org/) — in 0.x, a minor bump may be breaking.

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

[0.3.0]: https://github.com/vardister/Kraken.jl/releases/tag/v0.3.0
