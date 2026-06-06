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
        demes = instances(demeset.DemeSet)
        new{demeset.DemeSet,length(demes)}(
            Tuple(BitSet() for _ ∈ Base.OneTo(length(demes)))
        )
    end
    Coloring(y::Coloring{D,N}) where {D,N} = new{D,N}(copy(y.cols))
end

copy(c::NTuple{N,BitSet}) where N = Tuple(copy(i) for i in c)

copy(y::Coloring) = Coloring(y)

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
    @assert Int(b) ∈ y[i] "lineage $b is not in '$i' deme"
    delete!(y[i],Int(b))
    push!(y[j],Int(b))
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
    @assert Int(b) ∈ y[i] "lineage $b is not in '$i' deme"
    delete!(y[i],Int(b))
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
    @assert Int(b) ∈ y[i] "lineage $b is not in '$i' deme"
    chop!(y,i,Int(b))
    push!(y[j],Int(c))
    ell(y)
end

"""
    fork!(y, i, b, j, c)

Fork lineage `b` in deme `i` to lineages given in the vector `c` in demes given respectively by the vector `j`.
"""
fork!(
    y::Coloring{D},
    i::D, b::Integer,
    j::Union{NTuple{N,D},Vector{D}},
    c::Union{NTuple{N,I},Vector{I}},
) where {N,D,I<:Integer} = begin
    @assert Int(b) ∈ y[i] "lineage $b is not in '$i' deme"
    delete!(y[i],Int(b))
    for k ∈ eachindex(j)
        push!(y[j[k]],Int(c[k]))
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
    push!(y[i],Int(b))
    ell(y)
end
