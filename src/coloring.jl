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

ell(y::BitSet) = length(y)

"""
    swap!(y, i, j, b)

Swap lineage `b` from deme `i` to deme `j` in coloring `y`.
"""
swap!(
    y::Coloring{D},
    i::D, j::D,
    b::Integer,
) where D = begin
    @assert b ∈ y[i]
    delete!(y[i],b)
    push!(y[j],b)
    ell(y)
end

"""
    chop!(y, i, b)

Chop lineage `b` from deme `i`.
"""
chop!(
    y::Coloring{D},
    i::D,
    b::Integer,
) where D = begin
    @assert b ∈ y[i]
    delete!(y[i],b)
    ell(y)
end

"""
    chop!(y, i, b, j, c)

Chop lineage `b` from deme `i`; add lineage `c` in deme `j`.
"""
chop!(
    y::Coloring{D},
    i::D,
    b::Integer,
    j::D,
    c::Integer,
) where D = begin
    chop!(y,i,b)
    push!(y[j],c)
    ell(y)
end

"""
    fork!(y, i0, i1, i2, b0, b1, b2)

Fork lineage `b0` in deme `i0` to lineages `b1`, `b2` in demes `i1`, `i2`, respectively.
"""
fork!(
    y::Coloring{D},
    i0::D, i1::D, i2::D,
    b0::Integer, b1::Integer, b2::Integer,
) where D = begin
    @assert b0 ∈ y[i0]
    delete!(y[i0],b0)
    push!(y[i1],b1)
    push!(y[i2],b2)
    ell(y)
end

"""
    fork!(y, i, b, j, c)

Fork lineage `b` in deme `i` to lineages given in the vector `c` in demes given respectively by the vector `j`.
"""
fork!(
    y::Coloring{D},
    i::D, j::Vector{D},
    b::Integer, c::Vector{<:Integer},
) where D = begin
    @assert b ∈ y[i]
    delete!(y[i],b)
    for k ∈ eachindex(j)
        push!(y[j[k]],c[k])
    end
    ell(y)
end

"""
    plant!(y, i, b)

Plant lineage `b` in deme `i` within coloring `y`.
"""
plant!(
    y::Coloring{D},
    i::D,
    b::Integer,
) where D = begin
    push!(y[i],b)
    ell(y)
end

