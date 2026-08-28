using DrWatson
@quickactivate

using Test
using ForwardDiff
using BenchmarkTools
using FiniteDiff

include("../src/Kraken.jl")
using .Kraken

#%% Pekeris environment

function run_pekeris(freq, params; type="full")
    cw, cb, ρw, ρb, depth = params
    ssp, layers, sspHS = pekeris_env(; c0=cw, cb=cb, ρ0=ρw, ρb=ρb, depth=depth)
    env = UnderwaterEnv(ssp, layers, sspHS)
    solution = kraken_jl(env, freq)
    if type == "full"
        return solution
    elseif type == "kr"
        return solution.kr
    elseif type == "mode"
        return solution.modes
    end
end

function run_pekeris_fortran(freq, params; type="full")
    cw, cb, ρw, ρb, depth = params
    ssp, layers, sspHS = pekeris_env(; c0=cw, cb=cb, ρ0=ρw, ρb=ρb, depth=depth)
    env = EnvKRAKEN(ssp, layers, sspHS, [0, 100], 25.0)
    solution = kraken(env, freq; n_modes=5)
    if type == "full"
        return solution
    elseif type == "kr"
        return solution["kr_real"]
    elseif type == "mode"
        return solution["modes"]
    end
end

#%% Time both implementations
params = [1500.0, 1600.0, 1000.0, 1500.0, 100.0]

@btime run_pekeris(100.0, params)
@btime run_pekeris_fortran(100.0, params)

#%% Compare the results
sol = run_pekeris(100.0, params)
sol_fortran = run_pekeris_fortran(100.0, params)

#%% Time the derivatives
function dkr_dparams_jl(params)
    dkr_dparams = ForwardDiff.jacobian(x -> run_pekeris(100.0, x; type="kr"), params)
    return dkr_dparams
end

function dkr_dfreq_jl(freq, params)
    dkr_dparams = ForwardDiff.derivative(x -> run_pekeris(x, params; type="kr"), freq)
    return dkr_dparams
end

@btime dkr_dparams_jl(params)
@btime dkr_dfreq_jl(100.0, params)

function dkr_dparams_fortran(params)
    dkr_dparams = FiniteDiff.finite_difference_jacobian(x -> run_pekeris_fortran(100.0, x; type="kr"), params)
    return dkr_dparams
end

@btime dkr_dparams_fortran(params)

function calculate_kr(freq, params)
    cw, cb, ρw, ρb, depth = params
    ssp, layers, sspHS = pekeris_env(; c0=cw, cb=cb, ρ0=ρw, ρb=ρb, depth=depth)
    env = UnderwaterEnv(ssp, layers, sspHS)
    props = AcousticProblemProperties(env, freq)
    cache = AcousticProblemCache(env, props)
    krs = find_kr(env, props, cache)
    return krs
end
function dkr_dparams_jl_early(freq, params)
    dkr_dparams = ForwardDiff.jacobian(x -> calculate_kr(freq, x), params)
    return dkr_dparams
end

function dkr_dfreq_jl_early(freq, params)
    dkr_dparams = ForwardDiff.derivative(x -> calculate_kr(x, params), freq)
    return dkr_dparams
end

@btime dkr_dparams_jl_early(100.0, params)
@btime dkr_dfreq_jl_early(100.0, params)

function dpsi_dparams_jl(params)
    dkr_dparams = ForwardDiff.jacobian(x -> run_pekeris(100.0, x; type="mode"), params)
    return dkr_dparams
end

function dpsi_dparams_fortran(params)
    dkr_dparams = FiniteDiff.finite_difference_jacobian(
        x -> run_pekeris_fortran(100.0, x; type="mode"), params; absstep=1
    )
    return dkr_dparams
end
