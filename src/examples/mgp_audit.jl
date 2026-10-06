# Static audit and validation of an MGPModel.
# Reads the event table only; never calls a hazard closure.

export audit_model, validate_model, EventAudit, ModelAudit

"""
    EventAudit

One row of a model audit: static facts about one `Event`.
Deme names are resolved from the `from`/`into` indices.
"""
struct EventAudit
    name       :: Symbol
    type       :: EventType
    regular    :: Bool
    observed   :: Bool
    delta      :: Vector{Pair{Symbol,Int}}
    has_hazard :: Bool
    r          :: Vector{Int}
    from       :: Int
    from_deme  :: Union{Symbol,Nothing}
    into       :: Vector{Int}
    into_demes :: Vector{Symbol}
end

"""
    ModelAudit

Structural audit report for an `MGPModel`: its compartments, its
lineage-carrying demes, and one `EventAudit` per event.
Returned by `audit_model`. Plain printable value; `show` renders a table.
"""
struct ModelAudit
    name         :: Symbol
    compartments :: Vector{Symbol}
    demes        :: Vector{Symbol}
    events       :: Vector{EventAudit}
end

"""
    audit_model(model::MGPModel) -> ModelAudit

Inspect `model`'s Population IR and return a structured `ModelAudit`
exposing, for every event: its semantic type, regular/singular and
observed flags, state transition Δ_u, whether a hazard α_u is present,
production vector r_u, and source/target deme wiring W_u (`from`/`into`,
resolved to deme names as well as raw indices).

Never calls a hazard closure and never throws; see `validate_model` for checks.

`show(io, MIME("text/plain"), audit)` renders a human-readable table.
"""
function audit_model(model::MGPModel)
    demes = model.demes
    events = EventAudit[]
    for ev in model.events
        from_deme = (1 <= ev.from <= length(demes)) ? demes[ev.from] : nothing
        into_demes = Symbol[(1 <= i <= length(demes)) ? demes[i] : Symbol("<invalid:$i>")
                             for i in ev.into]
        push!(events, EventAudit(
            ev.name, ev.type, ev.regular, ev.observed,
            ev.Δ, ev.hazard isa Function, ev.r,
            ev.from, from_deme, ev.into, into_demes,
        ))
    end
    ModelAudit(model.name, model.compartments, demes, events)
end

function Base.show(io::IO, ::MIME"text/plain", ev::EventAudit)
    tag = ev.regular ? "regular" : "singular"
    obs = ev.observed ? "observed" : "background"
    println(io, "  - ", ev.name, "  [", ev.type, ", ", tag, ", ", obs, "]")
    delta_str = isempty(ev.delta) ? "(none)" : join(("$(p.first)=$(p.second >= 0 ? "+" : "")$(p.second)" for p in ev.delta), ", ")
    println(io, "      Δ (state jump):     ", delta_str)
    println(io, "      α (hazard):         ", ev.has_hazard ? "present (closure)" : "MISSING")
    println(io, "      r (production):     ", ev.r)
    from_str = ev.from_deme === nothing ? "(none)" : "$(ev.from_deme) [deme #$(ev.from)]"
    println(io, "      W.from (source):    ", from_str)
    into_str = isempty(ev.into) ? "(none)" :
        join(("$(d) [deme #$(i)]" for (d, i) in zip(ev.into_demes, ev.into)), ", ")
    println(io, "      W.into (targets):   ", into_str)
end

Base.show(io::IO, ev::EventAudit) = show(io, MIME("text/plain"), ev)

function Base.show(io::IO, ::MIME"text/plain", a::ModelAudit)
    println(io, "ModelAudit(", a.name, ")")
    println(io, "  compartments: ", join(a.compartments, ", "))
    println(io, "  demes:        ", join(a.demes, ", "))
    println(io, "  events (", length(a.events), "):")
    for ev in a.events
        show(io, MIME("text/plain"), ev)
    end
end

Base.show(io::IO, a::ModelAudit) = show(io, MIME("text/plain"), a)

"""
    validate_model(model::MGPModel) -> Vector{String}

Run structural/semantic checks against `model` and return a list of
human-readable issue strings; an empty vector means the model passed every
check. Never throws.

Covers: duplicate compartment/deme/event names, every deme being a
declared compartment, every event's `from`/`into` deme indices being valid
indices into `model.demes`, every event's `r` vector length matching
`length(model.demes)`, every event's hazard being a callable `Function`,
and every Δ entry referencing a declared compartment once.

Per `EventType`, when `r` and `from` are in range:
- `BIRTH`: `sum(r) == 2`.
- `MIGRATION`: one `into` deme, and `r` is 1 there and 0 elsewhere.
- `DEATH`: `r` is all 0.
- `SAMPLE`: `r[from]` is 0 or 1, and `r` is 0 elsewhere.
- `NEUTRAL`: `from == 0`, no `into`, and `r` is all 0.
- Every event: for each deme `d`, Δ of that compartment equals
  `r[d] - (d == from)`. One lineage is kept per individual.

[`simulate`](@ref) throws `ArgumentError` when this returns any issue.
"""
function validate_model(model::MGPModel)
    issues = String[]
    ndemes = length(model.demes)
    comp_set = Set(model.compartments)

    allunique(model.compartments) ||
        push!(issues, "duplicate compartment names: $(model.compartments)")
    allunique(model.demes) ||
        push!(issues, "duplicate deme names: $(model.demes)")
    for d in model.demes
        d in comp_set || push!(issues, "deme `$d` is not a declared compartment")
    end
    allunique(ev.name for ev in model.events) ||
        push!(issues, "duplicate event names among model.events")

    for ev in model.events
        tag = "event `$(ev.name)`"

        ev.hazard isa Function ||
            push!(issues, "$tag: hazard (α_u) is not a callable Function")

        length(ev.r) == ndemes ||
            push!(issues, "$tag: r (production vector) has length $(length(ev.r)), " *
                           "expected $(ndemes) (= number of demes)")

        (ev.from == 0 || (1 <= ev.from <= ndemes)) ||
            push!(issues, "$tag: from index $(ev.from) out of range 1:$(ndemes)")

        for i in ev.into
            (1 <= i <= ndemes) ||
                push!(issues, "$tag: into index $(i) out of range 1:$(ndemes)")
        end

        for p in ev.Δ
            p.first in comp_set ||
                push!(issues, "$tag: Δ references unknown compartment `$(p.first)`")
        end

        allunique(p.first for p in ev.Δ) ||
            push!(issues, "$tag: Δ names a compartment more than once")

        if ev.type in (BIRTH, MIGRATION, DEATH, SAMPLE) && ev.from == 0
            push!(issues, "$tag: $(ev.type) event has from=0 (no source deme)")
        end
        if ev.type == MIGRATION && isempty(ev.into)
            push!(issues, "$tag: MIGRATION event has no into deme")
        end

        ## The rules below index `r` by deme and by `from`.
        (length(ev.r) == ndemes && 0 <= ev.from <= ndemes) || continue
        if ev.type == BIRTH
            sum(ev.r) == 2 ||
                push!(issues, "$tag: BIRTH has $(sum(ev.r)) products (sum of r), expected 2")
        elseif ev.type == MIGRATION
            length(ev.into) == 1 ||
                push!(issues, "$tag: MIGRATION has $(length(ev.into)) into demes, expected 1")
            if length(ev.into) == 1 && 1 <= only(ev.into) <= ndemes
                ev.r == [Int(j == only(ev.into)) for j in 1:ndemes] ||
                    push!(issues, "$tag: MIGRATION r = $(ev.r) is not 1 at its into deme and 0 elsewhere")
            end
        elseif ev.type == DEATH
            all(iszero, ev.r) || push!(issues, "$tag: DEATH r = $(ev.r), expected all 0")
        elseif ev.type == SAMPLE && ev.from >= 1
            ev.r[ev.from] in (0, 1) && all(iszero, ev.r[j] for j in 1:ndemes if j != ev.from) ||
                push!(issues, "$tag: SAMPLE r = $(ev.r); r[from] must be 0 or 1 and the rest 0")
        elseif ev.type == NEUTRAL
            ev.from == 0 && isempty(ev.into) && all(iszero, ev.r) ||
                push!(issues, "$tag: NEUTRAL must have from = 0, no into, and r all 0")
        end
        for (d, deme) in enumerate(model.demes)
            implied = ev.r[d] - (d == ev.from ? 1 : 0)
            stated = sum((p.second for p in ev.Δ if p.first == deme); init = 0)
            stated == implied ||
                push!(issues, "$tag: Δ changes deme `$deme` by $stated, but r and from imply $implied")
        end
    end

    issues
end
