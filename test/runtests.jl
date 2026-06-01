using PhyloPOMP
using Test
using Crayons

h1 = crayon"bold blue"
h2 = s -> crayon"!bold light_yellow"("- "*s)

@testset verbose=true "PhyloPOMP.jl" begin

    include("parse.jl")
    include("newick.jl")
    include("fsmarkov.jl")
    include("guide.jl")
    include("rcateg.jl")
    include("seir_naive.jl")
    include("seir_guided.jl")

end
