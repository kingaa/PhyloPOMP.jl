import Base: getindex, copy

"""
        Coloring(d)

The `Coloring` struct holds information on the local coloring of a
genealogy.  Construct it with the call `Coloring(d)`, where `d` is
a demeset (see [@demes](@ref)).
"""
mutable struct Coloring{D <: Enum, N}
    cols::NTuple{N,BitSet}
    Coloring(demeset::Module) = begin
        demes = instances(demeset.T)
        new{demeset.T,length(demes)}(
            Tuple(BitSet() for _ ∈ Base.OneTo(length(demes)))
        )
    end
end

copy(y::Coloring) = deepcopy(y)

getindex(y::Coloring{D}, i::D) where D = getindex(y.cols, Int(i))

ell(y::Coloring) = length.(y.cols)

ell(y::Coloring{D}, i::D) where D = length(y.cols[Int(i)])
