default:
	julia --project -e 'import Pkg; Pkg.instantiate(); Pkg.build(); Pkg.precompile()'
	julia --project=@. make.jl

.PHONY: default
