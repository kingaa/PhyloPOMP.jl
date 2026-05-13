using PhyloPOMP
using Test

@testset verbose=true "PhyloPOMP.jl" begin

    include("parse.jl")
    include("newick.jl")
    include("fsmarkov.jl")
    include("guide.jl")
    include("rcateg.jl")

end
