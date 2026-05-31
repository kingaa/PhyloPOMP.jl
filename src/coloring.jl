"""
        Coloring(d)

The `Coloring` struct holds information on the local coloring of a
genealogy.  Construct it with the call `Coloring(d)`, where `d` is
a demeset (see [@demes](@ref)).
"""
mutable struct Coloring{D <: Enum}
    cols::Vector{BitSet}
    Coloring(demeset::Module) = begin
        demes = instances(demeset.T)
        new{demeset.T}(
            [BitSet() for _ ∈ Base.OneTo(length(demes))]
        )
    end
end

import Base: getindex
getindex(y::Coloring{D}, i::D) where D = getindex(y.cols, Int(i))
