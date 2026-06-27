using EnumX: @enumx

"""
    @demes Name deme1 deme2 ...

Creates a new demeset, `Name`, which contains an enumeration of the demes
`deme1`, `deme2`, ....  A demeset is a `Module`, which contains the
enumeration and the enumeration type.
"""
macro demes(name, first, rest...)
    esc(:(@enumx T=DemeSet $name $first=1 $(rest...)))
end

"""
    Unstructured

A demeset for the unstructured (1-deme) case.
"""
@demes Unstructured default

"""
    name2enum(demes)

Returns a function that maps strings to demes. It will also interpret a given
integer according to the integer representation of the enumeration.
"""
name2enum(demes::Type{D}) where {D <: Enum} = begin
    d = merge(
        Dict(lowercase(String(Symbol(i))) => i for i ∈ instances(demes)),
        Dict("$(Int(i))" => i for i ∈ instances(demes)),
    )
    function (name::AbstractString)
        get(d,lowercase(name),missing)
    end
end

"""
    enum2name(demeset)

Returns a `Dict` that maps demes to strings.
"""
enum2name(demes::Type{D}) where {D <: Enum} = begin
    Dict(i=>String(Symbol(i)) for i ∈ instances(demes))
end

macro isademeset(d)
    esc(:(@assert $d isa Module && isdefined($d,:DemeSet) "construct a demeset using `@demes`"))
end
