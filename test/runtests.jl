using PhyloPOMP
using Test
using Crayons

h1 = crayon"bold blue"
h2 = s -> crayon"!bold light_yellow"("- "*s)

@testset verbose=true "PhyloPOMP.jl" begin

    include("parse.jl")
    include("newick.jl")
    include("cblv.jl")
    include("fsmarkov.jl")
    include("guide.jl")
    include("rcateg.jl")
    include("seir_macro_equivalence.jl")
    include("seir_simulate.jl")
    include("seir_naive.jl")
    include("seir_soft.jl")
    include("seir_guided.jl")
    include("seir_hard.jl")
    include("mers_naive.jl")
    include("mers_simulate.jl")
    include("mers_soft.jl")
    include("mers_guided.jl")
    include("mers_hard.jl")

end
