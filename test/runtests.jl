using PhyloPOMP
using Test

@testset verbose=true "PhyloPOMP.jl" begin

    include("parse.jl")
    include("fsmarkov.jl")
    include("guide.jl")

end
