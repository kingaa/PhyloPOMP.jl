import Base: show

Base.show(
    io::IO,
    g::Union{Genealogy,GenealNode,Guide,GuideNode};
    kwargs...,
) = begin
    print(io,pretty_string(g;kwargs...))
end

pretty_string(g::Genealogy; sigdigits=4) = begin
    "<genealogy on ["*
        "$(round(g.t0,sigdigits=sigdigits)),"*
        "$(round(g.time,sigdigits=sigdigits))]:\n" *
        join(
            map(eachindex(g)) do i
                "  $i: " * pretty_string(g[i],sigdigits=sigdigits)
            end,
            '\n'
        ) * ">"
end

pretty_string(p::GenealNode; sigdigits = 4) = begin
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
    "<$(lowercase(String(Symbol(p.type)))) " *
        "lineage=$(p.lineage) " *
        "deme=$(p.deme) " *
        "time=$(round(p.slate,sigdigits=sigdigits)) " *
        parent_string *
        "children=[$children_string]" * ">"
end

pretty_string(g::Guide; sigdigits = 4) = begin
    nodes = map(eachindex(g)) do i
        "  $i: "*pretty_string(g[i],sigdigits=sigdigits)
    end
    "<guide:\n" * join(nodes,'\n') * ">"
end

pretty_string(n::GuideNode; sigdigits = 4) = begin
    t1 = round(n.tbeg,sigdigits=sigdigits)
    t2 = round(n.tend,sigdigits=sigdigits)
    chillins = map(eachindex(n.chillins)) do i
        ell = n.chillins[i]
        prob = round.(n.present[:,i],sigdigits=sigdigits)
        "$ell=>$prob"
    end
    targs = map(eachindex(n.alllins)) do i
        ell = n.alllins[i]
        prob = round.(n.target[:,i],sigdigits=sigdigits)
        "$ell=>$prob"
    end
    "<$(lowercase(String(Symbol(n.type)))): t ∈ [$t1,$t2] " *
        "parlin=$(n.parlin)" *
        " children:{" * join(chillins,',') *
        "} targets:{" * join(targs,',') * "}>"
end
