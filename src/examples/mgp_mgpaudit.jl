# mgpaudit / @mgpaudit: per-mark derivation report.
# Chain: audit_model -> full_transitions -> reduce_event_indicator
#        -> classify_filter_terms -> explain.
# Included after mgp_explain.jl because the report structs name its types.

export mgpaudit, @mgpaudit, default_audit_state,
       MarkAudit, MarkAuditInScope, MarkAuditOutOfScope, MgpAuditReport

"""
    default_audit_state(model::MGPModel) -> (ℓ::Vector{Int}, n::Vector{Int})

Default `(ℓ, n)` for `mgpaudit`:
- SEIR: `([2,2], [6,5])`
- MERS: `([2,2], [5,5])`
- other models: `(fill(2, D), fill(5, D))`, `D = length(model.demes)`

Large enough that each BIRTH mark able to fork reaches a Fork saturation.
"""
function default_audit_state(model::MGPModel)
    if model.name == :SEIR
        return (Int[2, 2], Int[6, 5])
    elseif model.name == :MERS
        return (Int[2, 2], Int[5, 5])
    else
        D = length(model.demes)
        return (fill(2, D), fill(5, D))
    end
end

"""
    MarkAudit

One event's row in an `MgpAuditReport`: `MarkAuditInScope` or
`MarkAuditOutOfScope`.
"""
abstract type MarkAudit end

"""
    MarkAuditInScope

Derivation of a BIRTH/MIGRATION event at `(ℓ, n)`: `event_audit`, the reduced
explanations, the `FilterSpec`, and one `TermExplanation` per term.
"""
struct MarkAuditInScope <: MarkAudit
    event_audit :: EventAudit
    ℓ           :: Vector{Int}
    n           :: Vector{Int}
    reduced     :: Vector{ReducedExplanation}
    spec        :: FilterSpec
    terms       :: Vector{TermExplanation}
end

"""
    MarkAuditOutOfScope

Record for a DEATH/SAMPLE/NEUTRAL event. `note` says why the event is outside
the saturation chain.
"""
struct MarkAuditOutOfScope <: MarkAudit
    event_audit :: EventAudit
    note        :: String
end

"""
    MgpAuditReport

Result of `mgpaudit`: `model_name`, the `(ℓ, n)` state, and one `MarkAudit` per
event in model order.
"""
struct MgpAuditReport
    model_name :: Symbol
    ℓ          :: Vector{Int}
    n          :: Vector{Int}
    marks      :: Vector{MarkAudit}
end

const OUT_OF_SCOPE_NOTE = (ev_name, ev_type) -> """
out of scope for the saturation chain: $(ev_type) events are decay/rate-driven, \
not saturation-driven.\
"""

"""
    mgpaudit(model::MGPModel; ℓ = nothing, n = nothing) -> MgpAuditReport

Derive every event of `model` at `(ℓ, n)` (`default_audit_state(model)` if
omitted). Returns an `MgpAuditReport`.
Throws `ArgumentError` unless `length(ℓ) == length(n) == length(model.demes)`
and `ℓ .<= n`.
"""
function mgpaudit(model::MGPModel; ℓ::Union{Nothing,AbstractVector{<:Integer}} = nothing,
                   n::Union{Nothing,AbstractVector{<:Integer}} = nothing)
    dℓ, dn = default_audit_state(model)
    ℓ = ℓ === nothing ? dℓ : collect(Int, ℓ)
    n = n === nothing ? dn : collect(Int, n)

    D = length(model.demes)
    length(ℓ) == D ||
        throw(ArgumentError("mgpaudit: ℓ has length $(length(ℓ)), expected $D " *
                             "(= length(model.demes))"))
    length(n) == D ||
        throw(ArgumentError("mgpaudit: n has length $(length(n)), expected $D " *
                             "(= length(model.demes))"))
    all(ℓ[d] <= n[d] for d in 1:D) ||
        throw(ArgumentError("mgpaudit: ℓ must be <= n elementwise (tracked lineages " *
                             "cannot exceed total population); got ℓ=$ℓ, n=$n"))

    audit = audit_model(model)
    marks = MarkAudit[]
    for (ev, ev_audit) in zip(model.events, audit.events)
        if ev.type in (BIRTH, MIGRATION)
            rts  = reduced_transitions(ev, ℓ, n)
            spec = classify_filter_terms(ev, rts)
            reduced_expl = [explain(rt) for rt in rts]
            all_terms = FilterTerm[spec.regular; spec.singular; spec.outflow_imbalance]
            term_expl = [explain(t) for t in all_terms]
            push!(marks, MarkAuditInScope(ev_audit, ℓ, n, reduced_expl, spec, term_expl))
        else
            push!(marks, MarkAuditOutOfScope(ev_audit, OUT_OF_SCOPE_NOTE(ev.name, ev.type)))
        end
    end

    MgpAuditReport(model.name, ℓ, n, marks)
end

function Base.show(io::IO, ::MIME"text/plain", m::MarkAuditInScope)
    ev = m.event_audit
    println(io, "  * ", ev.name, "  [", ev.type, ", r=", ev.r, "]  -- IN SCOPE" *
                 " (BIRTH/MIGRATION derivation)")
    println(io, "      state: ℓ=", m.ℓ, "  n=", m.n)
    println(io, "      reduced transitions: ", length(m.reduced))
    for r in m.reduced
        println(io, "        - ", r.kind, "  key=", r.key, "  Φ_u=", r.Φ,
                     "  (", length(r.members), " full transition(s) collapsed)")
    end
    nreg = count(t -> t.bucket == :regular_flow, m.terms)
    nsin = count(t -> t.bucket == :singular_flow, m.terms)
    noi  = count(t -> t.bucket == :outflow_imbalance, m.terms)
    println(io, "      classification: RegularFlow=", nreg,
                 "  SingularFlow=", nsin, "  OutflowImbalanceTerm=", noi)
    for t in m.terms
        println(io, "        - [", t.bucket, "] Φ=", t.Φ,
                     t.bucket == :outflow_imbalance ?
                         "  mechanism=$(t.mechanism) reason=$(t.reason) (NOT λ)" : "")
    end
end
Base.show(io::IO, m::MarkAuditInScope) = show(io, MIME("text/plain"), m)

function Base.show(io::IO, ::MIME"text/plain", m::MarkAuditOutOfScope)
    ev = m.event_audit
    println(io, "  * ", ev.name, "  [", ev.type, "]  -- OUT OF SCOPE")
    println(io, "      ", m.note)
end
Base.show(io::IO, m::MarkAuditOutOfScope) = show(io, MIME("text/plain"), m)

function Base.show(io::IO, ::MIME"text/plain", r::MgpAuditReport)
    println(io, "MgpAuditReport(", r.model_name, ")  state: ℓ=", r.ℓ, "  n=", r.n)
    println(io, "  marks (", length(r.marks), "):")
    for m in r.marks
        show(io, MIME("text/plain"), m)
    end
end
Base.show(io::IO, r::MgpAuditReport) = show(io, MIME("text/plain"), r)

"""
    @mgpaudit ModelName
    @mgpaudit ModelName ℓ=[...] n=[...]

Expands to `mgpaudit(ModelName; ℓ=..., n=...)` and returns its `MgpAuditReport`.
Further `key=value` arguments are forwarded as keyword arguments.
"""
macro mgpaudit(exprs...)
    isempty(exprs) &&
        throw(ArgumentError("@mgpaudit requires a model argument, e.g. `@mgpaudit SEIR`"))
    model_expr = exprs[1]
    kwexprs = exprs[2:end]
    for kw in kwexprs
        (kw isa Expr && kw.head === :(=)) ||
            throw(ArgumentError("@mgpaudit: expected `key=value` arguments after the " *
                                 "model, got `$(kw)`"))
    end
    return esc(:(mgpaudit($model_expr; $(kwexprs...))))
end
