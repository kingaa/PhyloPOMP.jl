# Event table of a model with the weights the generated filter uses, as Markdown or LaTeX.

export filter_table

"""
    filter_table(model; format = :markdown) -> String

One row per event of `model`: its name and type, its rate as written in `@mgp`, its change to the population,
what it does to the genealogy, and the weights of the generated filter (`mgp_filter_pomp`). Weights are written
in the counts `n` and tracked lineages `ℓ` of each deme; a prime marks a count after the event. Each formula
for a BIRTH or MIGRATION event is checked against the IR (`full_transitions`) on all states with up to 6 hosts
per deme before the table is returned; a mismatch throws.

`format = :latex` gives a `tabular` (Greek letters as `\\beta` etc.).
"""
function filter_table(model::MGPModel; format::Symbol = :markdown)
    format in (:markdown, :latex) || throw(ArgumentError("filter_table: format must be :markdown or :latex"))
    D = model.demes
    rows = Vector{NTuple{6,String}}()
    for ev in model.events
        _check_table_formulas(ev, length(D))
        rate = ev.rate === nothing ? "(closure)" : string(ev.rate)
        pop = isempty(ev.Δ) ? "none" : join(("$(s) $(k > 0 ? "+" : "")$(k)" for (s, k) in ev.Δ), ", ")
        tree, w = _table_row(ev, D)
        push!(rows, (string(ev.name), string(ev.type), rate, pop, tree, w))
    end
    format == :markdown ? _markdown(rows) : _latex(rows)
end

function _table_row(ev::Event, D)
    a = ev.from
    X = a > 0 ? string(D[a]) : ""
    nX(p = "") = "n$(p)_$X"; lX = "ℓ_$X"
    if ev.type == BIRTH
        d = _cross_dest(ev)
        if d == 0
            return "fork in $X",
                   "between nodes, no fork: 1 − C($lX,2)/C($(nX("′")),2); at a branch node: 1/C($(nX("′")),2)"
        end
        Y = string(D[d])
        return "fork $X → $X, $Y",
               "no lineage moves: 1 − ℓ_$Y/n′_$Y; one tracked lineage moves $X → $Y: (1 − (ℓ_$X − 1)/n_$X)/n′_$Y; " *
               "at a branch node: 1/(n′_$X n′_$Y)"
    elseif ev.type == MIGRATION
        d = _cross_dest(ev)
        Y = string(D[d])
        return "host moves $X → $Y",
               "untracked host (share 1 − ℓ_$X/n_$X): 1 − ℓ_$Y/n′_$Y; tracked lineage moves: 1/n′_$Y"
    elseif ev.type == DEATH
        return "lineage in $X ends", "only an untracked host; decay α·ℓ_$X/n_$X"
    elseif ev.type == SAMPLE
        destructive = any(p -> p.first == D[a] && p.second < 0, ev.Δ)
        node = destructive ? "at a sample tip: α" :
               "at a sample tip: α·(n_$X − ℓ′_$X)/n_$X; at a sampled ancestor: α/n_$X"
        return destructive ? "sample in $X, host removed" : "sample in $X, host stays", "$node; decay α"
    end
    "none", "no lineage involved"
end

## The BIRTH and MIGRATION formulas of `_table_row`, against the IR.
function _check_table_formulas(ev::Event, ndemes::Integer)
    ev.type in (BIRTH, MIGRATION) || return nothing
    a = ev.from; d = _cross_dest(ev)
    ## every state with up to 6 hosts in each deme
    for nv in Iterators.product(ntuple(_ -> 0:6, ndemes)...), lv in Iterators.product((0:k for k in nv)...)
        n = collect(nv); ℓ = collect(lv)
        all(n .>= ev.r) || continue          # n is the post-event count
        ok = true
        if ev.type == BIRTH && d == 0
            ok &= _no_move_target_ir(ev, ℓ, n) ≈ 1 - ℓ[a] * (ℓ[a] - 1) / (n[a] * (n[a] - 1))
            if ℓ[a] >= 2
                ts = full_transitions(ev, ℓ, n)
                ok &= Float64(_phi(ts, t -> t isa ForkTransition)) ≈ 2 / (n[a] * (n[a] - 1))
            end
        else
            ok &= _no_move_target_ir(ev, ℓ, n) ≈ 1 - ℓ[d] / n[d]
            if ℓ[d] >= 1                     # ℓ after a tracked lineage moved into d
                cr = ev.type == BIRTH ? (1 - ℓ[a] / n[a]) / n[d] : 1 / n[d]
                ok &= _cross_phi_ir(ev, d, ℓ, n) ≈ cr
            end
            if ev.type == BIRTH && ℓ[a] >= 1 && ℓ[d] >= 1
                ts = full_transitions(ev, ℓ, n)
                ok &= Float64(_phi(ts, t -> t isa ForkTransition)) ≈ 1 / (n[a] * n[d])
            end
        end
        ok || error("filter_table: a formula for event `$(ev.name)` does not match the IR at ℓ = $ℓ, n = $n")
    end
    nothing
end

function _markdown(rows)
    head = "| Event | Type | Rate | Population | Genealogy | Filter weights |\n|---|---|---|---|---|---|\n"
    head * join(("| " * join(r, " | ") * " |" for r in rows), "\n") * "\n"
end

const _GREEK = Dict('α' => "\\alpha", 'β' => "\\beta", 'γ' => "\\gamma", 'δ' => "\\delta", 'η' => "\\eta",
                    'κ' => "\\kappa", 'λ' => "\\lambda", 'μ' => "\\mu", 'ω' => "\\omega", 'ψ' => "\\psi",
                    'χ' => "\\chi", 'σ' => "\\sigma", 'ℓ' => "\\ell", 'θ' => "\\theta", 'ρ' => "\\rho",
                    '−' => "-", '′' => "'", '→' => "\\to")
function _tex(s::AbstractString)
    out = IOBuffer()
    for c in s
        if haskey(_GREEK, c)
            print(out, "\$", _GREEK[c], "\$")
        elseif c == '_'
            print(out, "\\_")
        else
            print(out, c)
        end
    end
    replace(String(take!(out)), "\$\$" => "")
end
_latex(rows) = "\\begin{tabular}{llp{3cm}p{3cm}p{3cm}p{6cm}}\nEvent & Type & Rate & Population & Genealogy & Filter weights \\\\\n\\hline\n" *
    join((join(_tex.(r), " & ") * " \\\\" for r in rows), "\n") * "\n\\end{tabular}\n"
