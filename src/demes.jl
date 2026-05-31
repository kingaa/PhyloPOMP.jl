"""
    @demes Name deme1 deme2 ...

Creates a new demeset, `Name`, which contains an enumeration of the demes `deme1`, `deme2`, ....  A demeset is a `Module`, which contains the enumeration and the enumeration type.
"""
macro demes(name, first, rest...)
    expr = :(@enumx $name $first=1 $(rest...))
    esc(:@eval $expr)
end

@demes Unstructured default

"""
    name2enum(demetype)

Returns a function that maps strings to demes.
"""
name2enum(demes::Type{D}) where {D <: Enum} = begin
    d = Dict(lowercase(String(Symbol(i)))=>i for i ∈ instances(demes))
    function (name::AbstractString)
        get(d,lowercase(name),missing)
    end
end

"""
    enum2name(demetype)

Returns a `Dict` that maps demes to strings.
"""
enum2name(::Val{D}) where {D <: Enum} = begin
    Dict(i=>String(Symbol(i)) for i ∈ instances(D))
end
