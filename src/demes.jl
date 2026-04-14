"""
    @demes Name deme1 deme2 ...

Creates a new type, `Name`, which is a (1-based) enumeration of the demes
`deme1`, `deme2`, ....
"""
macro demes(name, first, rest...)
    expr = :(@enum $name $first=1 $(rest...))
    esc(:@eval $expr)
end

@demes Unstructured default

name2enum(demes::Type{D}) where {D <: Enum} = begin
    d = Dict(lowercase(String(Symbol(i)))=>i for i ∈ instances(demes))
    function (name::AbstractString)
        get(d,lowercase(name),missing)
    end
end

enum2name(::Val{D}) where {D <: Enum} = begin
    Dict(i=>String(Symbol(i)) for i ∈ instances(D))
end
