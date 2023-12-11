using Parameters
using DSP
using ProgressMeter
using Interpolations
using QuadGK
using FFTW

libpath = @__DIR__

fillnan(x) = isnan(x) ? zero(x) : x

"""
`env = Env(ssp, sspHS, bottom, n_receivers, z_receivers)`

Creates an underwater environment in the form of a sound speed profile (`ssp`), a
half-space sound speed profile (`sspHS`), and a bottom profile (`bottom`).

The environment is used to compute the propagation of sound waves using KRAKEN.

Required arguments:
- `ssp`: sound speed profile
- `sspHS`: half-space sound speed profile
- `bottom`: bottom profile
- `n_receivers`: number of receivers
- `z_receivers`: depth of receivers


Optional keywords:
- `n_layers`: number of layers in the environment (default: 2)
- `note1`: note for the sound speed profile (default: "NVW")
- `note2`: note for the half-space sound speed profile (default: "A")
- `bsig`: bottom interfacial roughness (default: 0)
- `n_sources`: number of sources (default: 1)
- `z_sources`: depth of sources (default: 500)
- `n_bc`: number of boundary conditions (default: size(ssp, 1))
"""
@with_kw struct Env
    n_layers::Int = 2

    note1::String = "NVW"

    note2::String = "A"
    bsig::Int = 0 # bottom interfacial roughness

    ssp::Matrix{Float64}
    sspHS::Matrix{Float64}
    bottom::Matrix{Float64}

    n_sources::Int = 1
    z_sources::Float64 = 500

    n_receivers::Int
    z_receivers::Matrix{Float64}
    n_bc::Int = size(ssp, 1)
end

"""
    ssp, b, sspHS = env_builder(; hw = 71.0, cw = 1471.0, ρw = 1.0, αw = 0.0, h1 = 10.0, c1 = 1500.0,
        ρ1 = 1.6, α1 = 0.025, cb = 1900.0, ρb = 2.0, αb = 0.25, type = "1layer_constant")

Helper function to build an underwater environment for use with creating an Env object.

Optional keywords:
- `hw`: water depth (default: 71)
- `cw`: water sound speed (default: 1471)
- `ρw`: water density (default: 1)
- `αw`: water attenuation (default: 0)
- `h1`: sediment depth (default: 10)
- `c1`: sediment sound speed (default: 1500)
- `ρ1`: sediment density (default: 1.6)
- `α1`: sediment attenuation (default: 0.025)
- `cb`: bottom sound speed (default: 1900)
- `ρb`: bottom density (default: 2)
- `αb`: bottom attenuation (default: 0.25)
- `type`: type of environment to build (default: "1layer_constant")

    - "1layer_constant": 1 layer, constant sound speed
"""
function env_builder(; hw=71.0, cw=1471.0, ρw=1.0, αw=0.0, h1=10.0, c1=1500.0,
    ρ1=1.6, α1=0.025, cb=1900.0, ρb=2.0, αb=0.25, type="1layer_constant")
    d0 = hw
    d1 = hw + h1

    b0 = [0.0 0.0 d0]
    ssp0 = [0.0 cw 0.0 ρw αw 0.0
        d0 cw 0.0 ρw αw 0.0]

    b1 = [0.0 0.0 d1]
    ssp1 = [d0 c1 0.0 ρ1 α1 0.0
        d1 c1 0.0 ρ1 α1 0.0]

    ssp_bhs = [d1 cb 0.0 ρb αb 0.0]
    ssp_ths = [0.0 343.0 0.0 0.00121 0.0 0.0]

    ssp = [ssp0; ssp1]
    sspHS = [ssp_ths; ssp_bhs]
    b = [b0; b1]

    return ssp, b, sspHS
end

"""
    kraken(env::Env, freq=15.0; n_modes=5, range_max=5e3, c_low=0.0, c_high=nothing)

Call KRAKEN to compute the propagation of sound waves in an underwater environment.

Required arguments:
- `env`: underwater environment

Optional keywords:
- `freq`: frequency (default: 15)
- `n_modes`: number of modes (default: 5)
- `range_max`: maximum range (default: 5e3)
- `c_low`: minimum sound speed (default: 0)
- `c_high`: maximum sound speed (default: maximum of SSP and SSPHS)
"""
function kraken(env::Env,
    freq=15.0;
    n_modes=5,
    range_max=5e3,
    c_low=0.0,
    c_high=nothing)
    if c_high === nothing
        c_high = maximum([maximum(env.ssp[:, 2]), maximum(env.sspHS[:, 2])])
    end

    c_low_high = [c_low c_high]
    @unpack n_layers, note1, n_bc, note2, bsig, ssp, sspHS, bottom, n_sources, z_sources, n_receivers, z_receivers = env

    cg, cp, kr_real, kr_imag, zm, modes = call_kraken(n_modes,
        freq,
        n_layers,
        [UInt8(x) for x in env.note1],
        bottom,
        n_bc,
        ssp,
        [UInt8(x) for x in env.note2],
        bsig,
        sspHS,
        c_low_high,
        range_max,
        n_sources,
        z_sources,
        n_receivers,
        z_receivers)
    return Dict("cg" => cg,
        "cp" => cp,
        "kr_real" => kr_real,
        "kr_imag" => kr_imag,
        "zm" => zm,
        "modes" => modes)
end

"""
    call_kraken(nm, frq, nl, note1, b, nc, ssp, note2, bsig, sspHS, clh, rng, nsr, zsr, nrc, zrc)

Direct call to KRAKEN library.
!!! To be used with `kraken` function only. Not recommened to use as is.

Required arguments:
- `nm`: number of modes
- `frq`: frequency
- `nl`: number of layers
- `note1`: note for SSP
- `b`: bottom profile
- `nc`: number of boundary conditions
- `ssp`: SSP
- `note2`: note for SSPHS
- `bsig`: bottom interfacial roughness
- `sspHS`: SSPHS
- `clh`: minimum and maximum sound speed
- `rng`: maximum range
- `nsr`: number of sources
- `zsr`: depth of sources
- `nrc`: number of receivers
- `zrc`: depth of receivers

"""
function call_kraken(nm, frq, nl, note1, b, nc, ssp, note2, bsig, sspHS, clh, rng, nsr, zsr, nrc, zrc)
    nz = nsr + nrc

    cg = zeros(1, nm)
    cp = zeros(1, nm)
    kr_real = zeros(1, nm)
    kr_imag = zeros(1, nm)
    zm = zeros(nz, 1)
    modes = zeros(nz, nm)

    ccall((:kraken_, "$libpath/kraken.dylib"),
        Nothing,
        (Ref{Int}, Ref{Float64}, Ref{Int}, Ref{UInt8}, Ref{Float64}, Ref{Int}, Ref{Float64}, Ref{UInt8}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Int}, Ref{Float64}, Ref{Int}, Ref{Float64}, Ref{Int}, Ref{Float64},
            Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
        nm, frq, nl, note1, b, nc, ssp, note2, bsig, sspHS, clh, rng, nsr, zsr, nrc, zrc, nz, cg, cp, kr_real, kr_imag, zm, modes)
    return cg, cp, kr_real, kr_imag, zm, modes
end

function pf_signal(env, ranges, zs, zr; T=2, fs=1000, n_modes=41)
    freqs = range(1, fs / 2; step=1 / T)
    rfft_freqs = rfftfreq(T * fs, fs)

    rho0 = 1000  # density of water (kg/m³)
    cw = env.ssp[1, 2]
    pf_signal = zeros(ComplexF64, length(ranges), length(rfft_freqs))
    for freq in freqs
        res = kraken(env, freq; n_modes=n_modes)
        kr = res["kr_real"] + 1im * res["kr_imag"] |> vec
        ϕ = res["modes"]

        _, zs_ind = findmin(abs.(res["zm"] .- zs)) # zs_ind in a CartesianIndex
        _, zr_ind = findmin(abs.(res["zm"] .- zr))
        ϕ_zs = ϕ[zs_ind[1], :]  # take first index of the CartesianIndex
        ϕ_zr = ϕ[zr_ind[1], :]
        # println(size(ϕ_zs))
        for (rr, r) in enumerate(ranges)
            t0 = r / cw - 0.1  # align window correctly in time
            Q = 1im*exp(-1im*pi/4) / (rho0*sqrt(8π*r))
            pf = @. Q*ϕ_zs*ϕ_zr*exp(-im*kr*r)/sqrt(kr)
            pf = fillnan.(pf)
            pf = sum(pf)
            pf_shifted = pf * exp(2im * pi * freq * t0)
            ind = findfirst(rfft_freqs .== freq)
            pf_signal[rr, ind] = pf_shifted
        end
    end
    return pf_signal
end
"""
    pf_adiabatic(freq, envs, ranges, zs, zr; n_modes = 41)

Calculate the adiabatic pressure field for a given frequency.

Required arguments:
- `freq`: frequency (Hz)
- `envs`: list of underwater environments (Vector{Env})
- `ranges`: ranges of the environments (m)
- `zs`: source depth (m)
- `zr`: receiver depth (m)

"""
function pf_adiabatic(freq, envs, ranges, zs, zr; n_modes=41)
    krs = []
    modes = []
    zm = Vector{Float64}
    @showprogress for (i, env) in enumerate(envs)
        res = kraken(env, freq; n_modes=n_modes)
        kr = res["kr_real"] + 1im * res["kr_imag"]
        push!(krs, kr)
        push!(modes, res["modes"])

        if i == 1
            zm = res["zm"]
        end
    end

    ## Interpolate wavenumbers in range for each mode
    krs_m = vcat(krs...)  # make krs_m be length(ranges) x length(modes)
    krs_m_interp = [linear_interpolation(ranges, col) for col in eachcol(krs_m)]
    krs_m_integral = [quadgk(krs_m_interp[i], 0, ranges[end])[1]
                      for i in eachindex(krs_m_interp)]

    ## Compute pressure(ω)
    rho0 = 1000  # density of water (kg/m³)
    Q = im * exp(-1im * pi / 4) / (rho0 * sqrt(8pi * ranges[end]))

    _, zs_ind = findmin(abs.(zm .- zs))
    _, zr_ind = findmin(abs.(zm .- zr))

    ϕ_zs = modes[1][zs_ind]
    ϕ_zr = modes[end][zr_ind]

    pressure_f = Q * ϕ_zs .* ϕ_zr .* exp.(-im * krs_m_integral) ./ sqrt.(vec(krs[end])) .|>
                 fillnan |> sum
    return pressure_f
end

"""
    pf_adiabatic_signal(envs, ranges, zs, zr; T=2, fs=1000, n_modes=41)

    Calculate the pressure field at a given receiver depth `zr` from a source at `zs` 
    and using the adiabatic approximation for a range-dependent environment.

    Required arguments:
    - `envs`: list of underwater environments (Vector{Env})
    - `ranges`: ranges of the environments (m)
    - `zs`: source depth (m)
    - `zr`: receiver depth (m)

    Optional keywords:
    - `T`: duration of the signal (sec) (default: 2)
    - `fs`: sampling frequency (Hz) (default: 1000)
    - `n_modes`: number of modes (default: 41)
"""
function pf_adiabatic_signal(envs, ranges, zs, zr; T=2, fs=1000, n_modes=41)
    freqs = range(1, fs / 2; step=1 / T)
    rfft_freqs = rfftfreq(T * fs, fs)

    cw = envs[1].ssp[1, 2]
    t0 = ranges[end] / cw - 0.1  # align window correctly in time
    pf_signal = zeros(ComplexF64, length(rfft_freqs))
    for freq in freqs
        pf = pf_adiabatic(freq, envs, ranges, zs, zr; n_modes=n_modes)
        pf_shifted = pf * exp(2im * pi * freq * t0)
        ind = findfirst(rfft_freqs .== freq)
        pf_signal[ind] = pf_shifted
    end
    return pf_signal
end
