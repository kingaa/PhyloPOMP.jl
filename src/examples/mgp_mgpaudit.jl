# mgp_mgpaudit.jl
# =============================================================================
# M06 Part 2: mgpaudit / @mgpaudit -- end-to-end per-mark derivation report.
#
#   A SEPARATE file from `mgp_audit.jl` (M01's `audit_model`/`validate_model`)
#   even though this is `audit_model`'s natural conceptual extension, purely
#   for an include-ORDER reason, not a design one: `MarkAuditInScope` below
#   carries fields typed `Vector{ReducedExplanation}`/`FilterSpec`/
#   `Vector{TermExplanation}`, all defined by files (`mgp_reduce.jl`,
#   `mgp_filter_ir.jl`, `mgp_explain.jl`) included in `Examples.jl` AFTER
#   `mgp_audit.jl`. Julia requires a struct's field types to already exist
#   at struct-DEFINITION time (unlike ordinary function bodies, which only
#   need their referenced globals to exist by CALL time) -- so this file is
#   included near the end of `Examples.jl`'s chain, after `mgp_explain.jl`,
#   while `mgp_audit.jl` stays exactly where M01 put it.
#
#   `mgpaudit` is an ordinary function (per this project's "macros are
#   syntax, not architecture" rule -- restated in `mgp_audit.jl`'s header).
#   For every BIRTH/MIGRATION event of a model it walks the FULL compiler
#   chain built by M01-M05(+M06's Part-0 rename):
#
#       audit_model (M01, EventAudit)
#         -> enumerate_saturations / kli_binomial_ratio (M02, mgp_phi.jl --
#            invoked internally by full_transitions, not reimplemented here)
#         -> full_transitions (M03, mgp_transitions.jl)
#         -> reduce_event_indicator (M04, mgp_reduce.jl)
#         -> classify_filter_terms (M05/M06, mgp_filter_ir.jl --
#            RegularFlow / SingularFlow / OutflowImbalanceTerm)
#         -> explain (M06 Part 1, mgp_explain.jl -- provenance walk back down
#            to the source Event for every term)
#
#   at ONE caller-supplied concrete (ℓ, n) lineage/population state (this
#   machinery is genuinely state-dependent -- see `default_audit_state`'s
#   docstring for why there is no state-independent "the" table, and this
#   milestone's choice of illustrative default per model).
#
#   For DEATH/SAMPLE/NEUTRAL events, `mgpaudit` does NOT invent a parallel
#   derivation -- it reports, explicitly, that the event is out of scope for
#   this chain, reusing the exact reasoning `full_transitions` (M03) and
#   `ChopTransition` (M03) already established: compatibility for this whole
#   event-type class is decay/rate-driven, not saturation-driven
#   (`mers_filter_suite.tex` lines 511-539), regardless of whether `r_u`
#   happens to be zero.
#
#   `@mgpaudit` is a THIN macro wrapper: it contains no derivation logic of
#   its own, only argument-forwarding to `mgpaudit`.
#
# Explicitly OUT of scope (same boundary as M05/M06 Part 0/1):
#   - Computing KLI's actual `lambda` (Eq. 47/B2) or any proposal/`pi_u`
#     logic -- `mgpaudit` reports the Filter IR's classification and
#     provenance only, exactly as M05/M06 built it.
#   - mgp_filter.jl's stubs and the 8 hand-coded filter modules (untouched).
# =============================================================================

export mgpaudit, @mgpaudit, default_audit_state,
       MarkAudit, MarkAuditInScope, MarkAuditOutOfScope, MgpAuditReport

"""
    default_audit_state(model::MGPModel) -> (ℓ::Vector{Int}, n::Vector{Int})

Pick a single, illustrative, per-model default lineage-count/population-count
state for `mgpaudit` to derive against when the caller does not supply one.

There is no state-INDEPENDENT saturation/Φ table -- `enumerate_saturations`
and `kli_binomial_ratio` (M02) are functions of `(ℓ, n)`, evaluated fresh at
whatever state the filter happens to be in at runtime. `mgpaudit` therefore
needs *a* concrete state to report against, and this function's job is only
to pick one that is illustrative (small, hand-checkable, and -- where
possible -- large enough in every lineage-carrying deme that a `Fork`-
classified saturation is actually reachable for every BIRTH mark of the
model, so the audit report shows all three Filter IR buckets rather than
silently omitting `OutflowImbalanceTerm` for lack of a large enough `ℓ`).

Choices made (documented, not hidden):
- `:SEIR` (demes `[:E, :I]`): `ℓ=[2,2]`, `n=[6,5]` -- the EXACT instance
  already hand-verified for `infection` in `test/kli_reduce_test.jl` /
  `test/kli_filter_ir_test.jl` (`Φ(noop)=8/15, Φ(cross)=1/10, Φ(fork)=1/30`),
  reused here rather than inventing new numbers. `ℓ_E=2 ≥ 1` and `ℓ_I=2 ≥ 1`
  so `infection`'s `(1,1)` Fork saturation is reachable; `progression`'s
  `r=(0,1)` structurally cannot Fork at any state (M04's finding), so its
  report is expected to show only `RegularFlow` terms -- not a bug.
- `:MERS` (demes `[Camel, Human]`, i.e. `[I_c, I_h]`): `ℓ=[2,2]`, `n=[5,5]`
  -- chosen so `ℓ_C=2` and `ℓ_H=2` both reach the `r_d=2` bound needed for
  `transmission_cc`'s (TCC) and `transmission_hh`'s (THH) same-deme Fork,
  while also `≥1` in both demes for `transmission_hc`/`transmission_ch`'s
  (THC/TCH) cross-deme Fork -- so all four BIRTH marks illustrate their
  `OutflowImbalanceTerm` bucket in the default report.
- Any other model: `ℓ = fill(2, D)`, `n = fill(5, D)` (`D =
  length(model.demes)`) -- a generic small illustrative fallback, using the
  same "ℓ=2 reaches typical r_d≤2 bounds" reasoning as above.

Callers needing a different, more targeted, or larger state should call
`mgpaudit(model; ℓ=..., n=...)` directly rather than relying on this default.
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

Abstract supertype for one event's row in an `MgpAuditReport`: either
`MarkAuditInScope` (BIRTH/MIGRATION, full M01-M06 derivation) or
`MarkAuditOutOfScope` (DEATH/SAMPLE/NEUTRAL, explicit "out of scope, here's
why" note -- never silently omitted).
"""
abstract type MarkAudit end

"""
    MarkAuditInScope

Full per-mark derivation for a BIRTH/MIGRATION event at the report's
`(ℓ, n)` state: `event_audit` (M01's static `EventAudit`), the `(ℓ, n)` state
it was derived at, the M04 `ReducedTransition`s (each wrapped in a
`ReducedExplanation`, M06 Part 1), the M05/M06 `FilterSpec` classification
itself, and the `TermExplanation` (M06 Part 1) of every one of its
`RegularFlow`/`SingularFlow`/`OutflowImbalanceTerm` terms.
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

An explicit "out of scope, here's why" record for a DEATH/SAMPLE/NEUTRAL
event -- `mgpaudit` never silently drops these from its report. `note`
explains the scope boundary, reusing the same reasoning already established
by `full_transitions`'s `ArgumentError` message and `ChopTransition`'s
docstring (`mgp_transitions.jl`, M03).
"""
struct MarkAuditOutOfScope <: MarkAudit
    event_audit :: EventAudit
    note        :: String
end

"""
    MgpAuditReport

Top-level structured result of `mgpaudit(model)`: the model name, the
`(ℓ, n)` state the in-scope marks were derived at, and one `MarkAudit` per
event of `model`, IN MODEL ORDER (so every event is covered, not just the
ones with something interesting to say).
"""
struct MgpAuditReport
    model_name :: Symbol
    ℓ          :: Vector{Int}
    n          :: Vector{Int}
    marks      :: Vector{MarkAudit}
end

const OUT_OF_SCOPE_NOTE = (ev_name, ev_type) -> """
out of scope for the saturation/φ_u derivation chain (M03-M06): $(ev_type) \
events' compatibility is decay/rate-driven, not saturation-driven, \
regardless of whether r_u happens to be zero -- see `full_transitions`'s \
ArgumentError message and `ChopTransition`'s docstring in \
src/examples/mgp_transitions.jl (M03), and mers_filter_suite.tex lines \
511-539 (removal's I_d-1>=ell_d outflow condition and the sub-threshold \
decay term; sampling's split between an observed-leaf singular update and a \
lambda-contributing background rate). Working out this event's actual \
decay/driver contribution is explicitly deferred to the later \
driver/boost/decay milestone (M07 per the master project plan); mgpaudit \
does not invent a parallel derivation for it.\
"""

"""
    mgpaudit(model::MGPModel; ℓ = nothing, n = nothing) -> MgpAuditReport

Ordinary function (macros are syntax, not architecture -- see this file's
header) implementing the FULL per-mark derivation `@mgpaudit` prints. For
every event of `model`, in model order:

- BIRTH/MIGRATION: walks `audit_model` (M01) -> `full_transitions` (M03,
  which internally calls M02's `enumerate_saturations`/`kli_binomial_ratio`)
  -> `reduce_event_indicator` (M04) -> `classify_filter_terms` (M05/M06) ->
  `explain` (M06 Part 1) at the report's `(ℓ, n)` state (`default_audit_state`
  if the caller does not supply one -- see its docstring for the concrete
  values chosen per model and why), and records a `MarkAuditInScope`.
- DEATH/SAMPLE/NEUTRAL: records a `MarkAuditOutOfScope` with an explicit
  note (never silently omitted) -- see `OUT_OF_SCOPE_NOTE`.

`ℓ`/`n`, if supplied, must have length `length(model.demes)` and satisfy
`ℓ .<= n` elementwise (the same invariant `kli_binomial_ratio`/
`enumerate_saturations` already assume); this is checked here and raises
`ArgumentError` rather than silently producing degenerate saturation tables.

This function performs NO new KLI math -- it is purely a presentation/
composition layer over the already-tested M01-M06 pipeline, per M06's
project-mandate scope.
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
    println(io, "      reduced transitions (M04): ", length(m.reduced))
    for r in m.reduced
        println(io, "        - ", r.kind, "  key=", r.key, "  Φ_u=", r.Φ,
                     "  (", length(r.members), " full transition(s) collapsed)")
    end
    nreg = count(t -> t.bucket == :regular_flow, m.terms)
    nsin = count(t -> t.bucket == :singular_flow, m.terms)
    noi  = count(t -> t.bucket == :outflow_imbalance, m.terms)
    println(io, "      classification (M05/M06): RegularFlow=", nreg,
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

THIN macro wrapper -- per this project's "macros are syntax, not
architecture" rule, ALL derivation logic lives in the ordinary function
`mgpaudit`; this macro only forwards its arguments to it. `ModelName` is
expanded, unevaluated, as the first positional argument (so `@mgpaudit SEIR`
becomes `mgpaudit(SEIR)`); any further `key=value` arguments are forwarded
verbatim as keyword arguments (so `@mgpaudit SEIR ℓ=[2,2] n=[6,5]` becomes
`mgpaudit(SEIR; ℓ=[2,2], n=[6,5])`).

Returns the `MgpAuditReport` `mgpaudit` returns -- at the REPL this
auto-displays via `Base.show`; in code, `report = @mgpaudit SEIR` captures it
like any other function call.
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
