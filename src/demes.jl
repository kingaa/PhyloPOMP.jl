@enum Unstructured default

name2enum(demes::Type{D}) where {D <: Enum} = begin
    Dict(lowercase(String(Symbol(i)))=>i for i ∈ instances(demes))
end

enum2name(::Val{D}) where {D <: Enum} = begin
    Dict(i=>String(Symbol(i)) for i ∈ instances(D))
end
