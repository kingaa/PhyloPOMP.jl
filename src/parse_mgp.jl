using Base.Iterators: drop
using .Meta: isexpr

parse_mgp_spec1(arg) = arg

parse_mgp_spec1(arg::Expr) = begin
    if isexpr(arg,:tuple) && all(isexpr.(arg.args,:(=)))
        dict = MgpDict()
        parse_mgp_spec2!(dict, arg.args...)
        dict
    else
        arg
    end
end

parse_mgp_spec2!(dict::MgpDict, args::Expr...) = begin
    for arg in args
        @assert isexpr(arg,:(=)) && arg.args[1] isa Symbol "misspecified argument $arg."
        dict[arg.args[1]] = parse_mgp_spec1(arg.args[2])
    end
    nothing
end

parse_mgp_spec!(dict::MgpDict, args::Expr...) = begin
    parse_mgp_spec2!(dict, args...)
    repair_mgp_spec!(dict)
    nothing
end

repair_mgp_spec!(dict) = begin
    @assert haskey(dict,:state) "state not specified"
    @assert haskey(dict,:enns) "enns not specified"
    @assert haskey(dict,:ells) "ells not specified"
    dict[:state] = make_named(dict[:state])
    nothing
end

peelblock(arg::Expr) = begin
    if isexpr(arg,:block,2) && arg.args[1] isa LineNumberNode
        arg.args[2]
    else
        arg
    end
end

make_named(arg::Expr) = begin
    @assert isexpr(arg,:tuple)
    :((;$(arg.args...),))
end

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
