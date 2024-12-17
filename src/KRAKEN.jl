module KRAKEN
	include("kraken_src.jl")

	export kraken
	export Env
	export EnvKRAKEN
	export env_builder
	export pf_adiabatic
	export pf_adiabatic_signal
	export pf_signal
	
	include("sus_signature.jl")
	export wilson_charge_jb

	# Plotting functions
	# include("kraken_plots.jl")
	# export plot_modes



	### New Julia implementation
	include("kraken_core.jl")

	export SampledSSP, SampledDensity
	export UnderwaterEnv, AcousticProblemProperties
	export prepare_vectors, bisection, solve_for_kr, inverse_iteration, det_sturm, kraken_jl, find_kr

	include("kraken_pekeris.jl")
	export PekerisUnderwaterEnv
end