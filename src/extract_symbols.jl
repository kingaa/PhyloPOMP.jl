using Base.Iterators: drop

extract_symbols(ex) = extract_symbols!(ex)

extract_symbols!(ex, syms = Set{Symbol}()) = begin
    if ex isa Symbol
        push!(syms, ex)
    elseif ex isa Expr
        nex = copy(ex)
        if ex.head ∈ (:call,:macrocall)
            nex.args = collect(drop(ex.args,1))
        end
        for arg in nex.args
            extract_symbols!(arg, syms)
        end
    end
    syms
end
