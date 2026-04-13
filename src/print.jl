import Base: show

Base.show(io::IO, g::Genealogy) = begin
    print(io,pretty_string(g))
end

Base.show(io::IO, p::GenealNode) = begin
    print(io,pretty_string(p))
end

pretty_string(g::Genealogy) = begin
    strings = map(i->"  "*pretty_string(g[i]),eachindex(g))
    "<genealogy on "*
        "[$(round(g.t0,sigdigits=6)),$(round(g.time,sigdigits=6))]:\n" *
        join(strings,'\n') * ">"
end

pretty_string(p::GenealNode) = begin
    parent_string = if isnothing(p.parent)
        ""
    else
        "parent=$(p.parent) "
    end
    children_string = if isempty(p.children)
        ""
    else
        join(map(i->"$i",p.children),',')
    end
    "<$(lowercase(String(Symbol(p.type)))) $(p.name): " *
        "lineage=$(p.lineage) " *
        "time=$(round(p.slate,sigdigits=6)) " *
        parent_string *
        "children=[$children_string]" * ">"
end
