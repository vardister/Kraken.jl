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
	include("kraken_plots.jl")
	export plot_modes



	### New Julia implementation
	include("src_standard_environments.jl")
	include("kraken_core.jl")


	export one_layer_test_dict_KRAKEN, one_layer_test_dict, pekeris_test_dict
	export SampledSSP, SampledDensity
	export UnderwaterEnv, AcousticProblemProperties
	export prepare_vectors, bisection, solve_for_kr, inverse_iteration, det_sturm, kraken_jl
end