const MgpDict = Dict{Symbol,Any}

parse_mgp_spec1(arg) = arg

parse_mgp_spec1(arg::Expr) = begin
    if arg.head == :tuple &&
        all(map(x -> (isa(x,Expr) && x.head==:(=)), arg.args))
        dict = MgpDict()
        parse_mgp_spec!(dict, arg.args...)
        dict
    else
        arg
    end
end

parse_mgp_spec!(dict::MgpDict, args::Expr...) = begin
    for arg in args
        @assert arg.head == :(=) && arg.args[1] isa Symbol "misspecified argument $arg."
        dict[arg.args[1]] = parse_mgp_spec1(arg.args[2])
    end
    nothing
end

make_named(arg::Expr) = begin
    @assert arg.head == :tuple
    :((;$(arg.args...),))
end

repair_mgp_spec!(dict) = begin
    @assert haskey(dict,:state) "state not specified"
    @assert haskey(dict,:enns) "enns not specified"
    @assert haskey(dict,:ells) "ells not specified"
    dict[:state] = make_named(dict[:state])
    nothing
end

peelblock(arg::Expr) = begin
    if arg.head == :block && length(arg.args)==2 &&
        arg.args[1] isa LineNumberNode
        arg.args[2]
    else
        arg
    end
end
