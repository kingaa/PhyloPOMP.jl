const _MOVE_TYPES = (:fork, :swap, :chop, :sample, :sample_remove, :none)

function _symbols(ex, label)
    xs = ex isa Expr && ex.head === :tuple ? ex.args : Any[ex]
    all(x -> x isa Symbol, xs) || error("@mgp: `$label` must be a tuple of names")
    Symbol[xs...]
end

function _rewrite(ex, comps, pars)
    if ex isa Symbol
        ex in comps && return :(x.$ex)
        ex in pars && return :(θ.$ex)
        return ex
    elseif ex isa Expr
        return Expr(ex.head, map(arg -> _rewrite(arg, comps, pars), ex.args)...)
    end
    ex
end

function _pairs(ex)
    terms = if ex isa Expr && ex.head === :tuple
        ex.args
    elseif ex isa Expr && ex.head === :(=)
        Any[ex]
    else
        error("@event: `pop` must contain changes such as `(S=-1, E=+1)`")
    end
    delta = Pair{Symbol,Int}[]
    for term in terms
        term isa Expr && term.head === :(=) ||
            error("@event: every population change must have the form `name=integer`")
        name, amount = term.args
        name isa Symbol && amount isa Integer ||
            error("@event: every population change must have the form `name=integer`")
        push!(delta, name => Int(amount))
    end
    delta
end

function _deme_index(deme, demes)
    deme isa Symbol || error("@event: a deme name must be a symbol")
    i = findfirst(==(deme), demes)
    isnothing(i) && error("@event: `$deme` is not one of the declared demes")
    i
end

function _arrow(ex)
    ex isa Expr && ex.head === :call && ex.args[1] === :(=>) ||
        error("@event: expected a move of the form `source => target`")
    ex.args[2], ex.args[3]
end

function _targets(first_target, rest)
    targets = first_target isa Expr && first_target.head === :tuple ?
        Any[first_target.args...] : Any[first_target]
    append!(targets, rest)
    all(x -> x isa Symbol, targets) || error("@event: move targets must be deme names")
    Symbol[targets...]
end

function _production(targets, demes)
    r = zeros(Int, length(demes))
    for target in targets
        r[_deme_index(target, demes)] += 1
    end
    r
end

function _decode_move(ex, demes)
    if ex === :none || (ex isa Expr && ex.head === :call && ex.args[1] === :none)
        return NEUTRAL, 0, Int[], zeros(Int, length(demes))
    end
    ex isa Expr && ex.head === :call ||
        error("@event: `move` must be one of $(_MOVE_TYPES)")
    op = ex.args[1]
    op in _MOVE_TYPES || error("@event: unknown move `$op`")

    if op === :fork
        source, first_target = _arrow(ex.args[2])
        targets = _targets(first_target, ex.args[3:end])
        from = _deme_index(source, demes)
        into_demes = copy(targets)
        parent = findfirst(==(source), into_demes)
        isnothing(parent) || deleteat!(into_demes, parent)
        into = unique(_deme_index(d, demes) for d in into_demes)
        return BIRTH, from, collect(into), _production(targets, demes)
    elseif op === :swap
        length(ex.args) == 2 || error("@event: `swap` takes one source-to-target move")
        source, target = _arrow(ex.args[2])
        targets = _targets(target, Any[])
        length(targets) == 1 || error("@event: `swap` must have one target deme")
        return MIGRATION, _deme_index(source, demes),
            [_deme_index(only(targets), demes)], _production(targets, demes)
    elseif op === :chop
        length(ex.args) == 2 || error("@event: `chop` takes one deme")
        return DEATH, _deme_index(ex.args[2], demes), Int[], zeros(Int, length(demes))
    elseif op === :sample
        length(ex.args) == 2 || error("@event: `sample` takes one deme")
        deme = ex.args[2]
        return SAMPLE, _deme_index(deme, demes), Int[], _production([deme], demes)
    elseif op === :sample_remove
        length(ex.args) == 2 || error("@event: `sample_remove` takes one deme")
        return SAMPLE, _deme_index(ex.args[2], demes), Int[], zeros(Int, length(demes))
    end
    error("@event: use `move=none`, not `none(...)`")
end

function _event_expr(args, comps, pars, demes)
    isempty(args) && error("@event: missing event name")
    name = first(args)
    name isa Symbol || error("@event: the event name must be a symbol")
    kv = Dict{Symbol,Any}()
    for assignment in args[2:end]
        assignment isa Expr && assignment.head === :(=) ||
            error("@event $name: expected `key=value`, got `$assignment`")
        key, value = assignment.args
        key isa Symbol || error("@event $name: invalid key `$key`")
        haskey(kv, key) && error("@event $name: duplicate key `$key`")
        kv[key] = value
    end
    haskey(kv, :rate) || error("@event $name: missing `rate`")
    haskey(kv, :move) || error("@event $name: missing `move`")
    allowed = Set((:rate, :pop, :move, :kind))
    unknown = setdiff(Set(keys(kv)), allowed)
    isempty(unknown) || error("@event $name: unknown keys $(collect(unknown))")

    kind = get(kv, :kind, :regular)
    kind in (:regular, :singular) ||
        error("@event $name: `kind` must be `regular` or `singular`")
    regular = kind === :regular
    event_type, from, into, r = _decode_move(kv[:move], demes)
    delta = _pairs(get(kv, :pop, Expr(:tuple)))
    hazard = :((x, θ) -> $(_rewrite(kv[:rate], comps, pars)))
    :(Event($(QuoteNode(name)), $delta, $hazard, $r, $event_type,
            $from, $into, $regular, $(!regular)))
end

"A marker parsed by `@mgp`; it is invalid outside an `@mgp` block."
macro event(args...)
    error("@event is only valid inside an @mgp block")
end

"""
    @mgp Name begin ... end

Turn a declarative Markov genealogy process specification into an auditable
`MGPModel` event table. The macro performs no likelihood calculation.
"""
macro mgp(name, block)
    name isa Symbol || error("@mgp: the model name must be a symbol")
    block isa Expr && block.head === :block || error("@mgp: expected a `begin ... end` block")

    declarations = Dict{Symbol,Vector{Symbol}}()
    event_calls = Expr[]
    for stmt in block.args
        stmt isa LineNumberNode && continue
        if stmt isa Expr && stmt.head === :(=)
            key, value = stmt.args
            key in (:compartments, :demes, :params) ||
                error("@mgp: unknown declaration `$key`")
            haskey(declarations, key) && error("@mgp: duplicate declaration `$key`")
            declarations[key] = _symbols(value, key)
        elseif stmt isa Expr && stmt.head === :macrocall && stmt.args[1] === Symbol("@event")
            push!(event_calls, stmt)
        else
            error("@mgp: unrecognized line `$stmt`")
        end
    end
    for key in (:compartments, :demes, :params)
        haskey(declarations, key) || error("@mgp: missing `$key` declaration")
    end

    comps = declarations[:compartments]
    demes = declarations[:demes]
    pars = declarations[:params]
    all(d -> d in comps, demes) || error("@mgp: every deme must also be a compartment")
    isempty(intersect(Set(comps), Set(pars))) ||
        error("@mgp: compartment and parameter names must be disjoint")
    events = [_event_expr(Any[call.args[3:end]...], comps, pars, demes) for call in event_calls]
    demes_name = Symbol(name, :Demes)

    esc(quote
        @demes $(demes_name) $(demes...)
        const $(name) = MGPModel($(QuoteNode(name)), $comps, $demes, [$(events...)])
    end)
end

@mgp SEIR begin
    compartments = (S, E, I, R)
    demes = (E, I)
    params = (β, σ, γ, ω, ψ, χ, N)

    @event infection   rate=β*S*I/N pop=(S=-1, E=+1) move=fork(I => E, I) kind=regular
    @event progression rate=σ*E     pop=(E=-1, I=+1) move=swap(E => I)    kind=regular
    @event recovery    rate=γ*I     pop=(I=-1, R=+1) move=chop(I)         kind=regular
    @event waning      rate=ω*R     pop=(R=-1, S=+1) move=none            kind=regular
    @event sampling    rate=ψ*I     pop=()           move=sample(I)      kind=singular
end
include("mgp_mers.jl")
