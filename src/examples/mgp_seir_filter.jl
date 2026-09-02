# mgp_seir_filter.jl
# =============================================================================
# M08: End-to-End SEIR Compiler — the first major project gate.
#
#   Assembles M01-M07's pieces (Event/MGPModel, full_transitions/
#   reduce_event_indicator (Phi_u), total_decay (lambda), driver/boost) into
#   an EXECUTABLE filter for SEIR specifically, built to be directly
#   comparable against `seir_naive.jl` (read-only oracle, never modified).
#
#   SCOPING DECISION (documented, not hidden): this file compiles the
#   REGULAR part of the filter (`compiled_regular_part!`, mirroring
#   `seir_naive.jl`'s `event_rates!`/`regular_part!`) via the compiler's own
#   machinery. The SINGULAR part (root-planting / sample-chop / branch-point
#   fork, `seir_naive.jl`'s `singular_part!`) is REUSED VERBATIM from
#   `NaiveSEIR.singular_part!` rather than re-derived, because M02-M06
#   already established (Gates 2/3, cited in the M08 handoff) that
#   `full_transitions`/`reduce_event_indicator` reproduce the same Phi_u
#   numbers `singular_part!` implicitly relies on -- re-deriving the singular
#   branch-point/chop/root logic from the IR here would duplicate already-
#   verified work without adding new coverage, and risks introducing NEW bugs
#   in a part of the pipeline this milestone was not asked to re-derive (its
#   job, per the M08 task, is the decay/driver/boost RECONCILIATION below).
#   `mgp_filter.jl`'s generic stubs (`kli_select`, `kli_decay`, `apply_move!`,
#   `singular_update!`) are intentionally left untouched (still `error()`):
#   filling them generically for every EventType would require committing to
#   and testing a canonical dispatch signature well beyond SEIR, which is out
#   of this milestone's scope (see the M08 task's own framing: "you do NOT
#   need to build a fully general compile_filter pass").
#
# =============================================================================
# THE DECAY DISCREPANCY (M07's flagged, unresolved question) — RESOLVED HERE
# =============================================================================
#
# M07 found `mgp_decay.jl`'s tex-literal DEATH formula
#     lambda_tex(DEATH) = alpha_u(x,theta) * 1{n_d <= ell_d}
# disagrees with `seir_naive.jl:120`/`mers_naive.jl:159-161`'s actual running
# `decay` accumulator, which adds `gamma_d*ell_d` UNCONDITIONALLY (every
# step, not just at the n_d==ell_d boundary). M07 flagged this, implemented
# the tex-literal formula anyway (matching the tex and `mgp_filter.jl`'s own
# `kli_decay` docstring), and left the reconciliation to M08.
#
# ISOLATED NUMERICAL EXPERIMENT (per the M08 task's mandate: resolve by a
# concrete, deterministic numerical comparison, not symbolic argument alone).
# `NaiveSEIR.event_rates!` was called DIRECTLY (not through the full
# stochastic simulator) at a fixed coloring (ellI=2 tracked I-lineages) and
# beta=sigma=omega=psi=chi=0 (gamma=1 only, isolating the DEATH/decay
# machinery), sweeping I from ellI up to ellI+5, and separately sweeping
# ellI itself from 0 to 4 (see this milestone's scratch probe,
# `decay_probe.jl`/`decay_probe2.jl` in the session scratchpad; every case
# below was additionally re-run and confirmed by hand):
#
#   I=5, ellI=2 (STRICTLY above threshold): alpha[5] (naive's REDUCED regular
#     proposal rate for recovery) = gamma*(I-ellI) = 3.0; naive's actual
#     `decay` return value = 2.0 = gamma*ellI (UNCONDITIONAL).
#     `mgp_decay.jl`'s tex-literal lambda at this state = 0.0 (I>ellI, so the
#     `1{I<=ellI}` gate is false) -- CONFIRMS the discrepancy exists exactly
#     as M07 described, using the real code, not a hand-transcription.
#   I=2=ellI=2 (boundary): alpha[5]=0.0; naive `decay`=2.0=gamma*ellI; tex
#     lambda=2.0 (gate true here) -- MATCHES at the boundary, exactly as
#     M07's own algebra predicted.
#
# RESOLUTION (first-principles derivation, then cross-checked against BOTH
# the tex's own text and `seir_naive.jl`'s numbers programmatically over many
# (I,ellI) combinations -- see `decay_probe2.jl`, "ALL MATCH: true"):
#
# `mgp_decay.jl`'s tex-literal lambda formula is CORRECT and needs NO
# change -- it is an exact transcription of `mers_filter_suite.tex`'s boxed
# "Decay lambda" equation (lines 833-836) and is the right formula for a
# DEATH event's contribution to the PDE-level lambda(t,x,y) alone.
#
# `seir_naive.jl`'s raw `decay` accumulator is NOT literally lambda(t,x,y).
# It is lambda PLUS an additional term that only exists because
# `seir_naive.jl`'s single-particle importance sampler PROPOSES regular
# DEATH jumps at a REDUCED rate (`alpha[5] = gamma*(I-ellI)`, i.e. "an
# untracked individual recovers" -- exactly `NaiveProposal`'s choice, M07
# Part B / M00's taxonomy), whereas `mers_filter_suite.tex`'s own "Assembled
# regular filter" section (lines 800-814) and its "reduced regular exit rate"
# notebox (lines 846-858) prescribe the REGULAR DEATH outflow/exit rate as
# the FULL, threshold-GATED hazard `gamma_d*I_d*1{I_d>ell_d}` (quoted
# directly: "above-threshold deaths: gamma_C I_C ind_{I_C>ell_C} + ...").
# These are two DIFFERENT rates for the same physical jump type:
#   alpha_u^reduced (naive's actual proposal rate)   = gamma_d*(n_d-ell_d)
#   alpha_u^tex-reg  (tex's prescribed regular exit)  = gamma_d*n_d*1{n_d>ell_d}
# Both are legitimate; a "naive"-style proposal is FREE to use a smaller,
# reduced rate for efficiency (it only ever proposes jumps guaranteed
# compatible with the data), but doing so means the DIFFERENCE between what
# the exact model's regular exit rate ought to be (tex's alpha_u^tex-reg) and
# what was actually proposed (alpha_u^reduced) is un-priced mass that MUST be
# charged somewhere -- and it is charged to `decay`, continuously, as a
# "leftover":
#     leftover_u(x,theta,ell,n) = alpha_u^tex-reg(x,theta,ell,n)
#                                - driver(alpha_u(x,theta), pi_u^naive(ell,n))
# where `pi_u^naive(ell,n) = (n_d-ell_d)/n_d` is exactly `NaiveProposal`'s
# selection fraction for this DEATH event (matching `alpha[5]`'s own
# `(I-ellI)/I` shape) and `driver` is M07's `driver(alpha,pi)=alpha*pi`
# (`mgp_proposal.jl`).
#
# ALGEBRAIC IDENTITY (why `naive_decay = lambda_tex + leftover` always, no
# case-split needed once expanded): writing `a=n_d-ell_d` (untracked count,
# `a>=0` by the model invariant),
#     lambda_tex   = gamma_d*n_d*1{a==0}
#     leftover     = gamma_d*n_d*1{a>0} - gamma_d*a
#     sum          = gamma_d*n_d*(1{a==0}+1{a>0}) - gamma_d*a
#                  = gamma_d*n_d - gamma_d*a  =  gamma_d*(n_d-a) = gamma_d*ell_d
# -- exactly `seir_naive.jl`/`mers_naive.jl`'s actual, unconditional
# `gamma_d*ell_d` decay contribution, for EVERY `n_d>=ell_d`, confirmed
# numerically above across (I,ellI) in [0,4]x[ellI,ellI+3] and MERS's own
# M07 instance (13/5 tex + 1 leftover(removal_c, I_C=5>ell_C=2) + 0
# leftover(removal_h, I_H=3=ell_H=3, boundary) = 18/5, M07's own cited naive
# total).
#
# CONCLUSION: `mgp_decay.jl` is UNCHANGED (still correct, as pure lambda).
# `compiled_decay` below adds the DEATH-specific `leftover` term (new to
# this milestone) so the ASSEMBLED filter's decay matches
# `seir_naive.jl`'s actual `decay` bookkeeping variable exactly -- this is
# the concrete, resolved answer to M07's flagged question.
# =============================================================================

export compiled_decay, compiled_event_rates!, compiled_regular_part!,
       compiled_filter_pomp

"""
    compiled_decay(model::MGPModel, x, θ, ℓ::AbstractVector{<:Integer},
                    n::AbstractVector{<:Integer}) -> Float64

`lambda(t,x,y)` (M07's `total_decay`, unmodified) PLUS the DEATH-event
"regular-proposal leftover" derived in this file's header, reproducing
`seir_naive.jl`'s/`mers_naive.jl`'s actual running `decay` bookkeeping
variable exactly (not just the tex's pure lambda) for a `NaiveProposal`-style
filter. `x`/`θ` are the population state/parameters; `ℓ`/`n` are per-deme
tracked/total counts, POST-event, in `model.demes` order (M02-M07 convention).
"""
function compiled_decay(model::MGPModel, x, θ,
                         ℓ::AbstractVector{<:Integer},
                         n::AbstractVector{<:Integer})
    total = Float64(total_decay(model, x, θ, ℓ, n))
    for event in model.events
        event.type == DEATH || continue
        d = event.from
        αfull = Float64(event.hazard(x, θ))
        nd, ℓd = n[d], ℓ[d]
        above = nd > ℓd
        αreduced = above ? αfull * (nd - ℓd) / nd : 0.0
        leftover = (above ? αfull : 0.0) - αreduced
        total += leftover
    end
    total
end

"""
    compiled_event_rates!(alpha, pi, cols, S, E, I, R; β, σ, γ, ω, ψ, χ, pop, model, _...) -> decay

Drop-in structural replacement for `NaiveSEIR.event_rates!`
(`seir_naive.jl:103-121`), same 6-slot `alpha`/`pi` layout and semantics
(1,2 = infection identity/cross; 3,4 = progression identity/cross; 5 =
recovery; 6 = waning), but every VALUE is sourced from the compiler:
`event.hazard` for `alpha`, `NaiveProposal`'s `(n-ell)/n` / `ell/n` split
(matching `seir_naive.jl`'s own `pi[1..4]` shapes, cited inline) for `pi`,
and `compiled_decay` (this file) for the return value. `χ` is accepted for
signature parity with `NaiveSEIR.filter_pomp`'s kwargs but -- per M07's
documented, confirmed representation gap (`mgp_decay.jl` header, "ONE
GENUINE GAP found") -- has no corresponding `Event` in `PhyloPOMP.SEIR`, so
it never enters `alpha`/`pi`/decay here (matches `total_decay`, which never
sees a `chi`-rated Event either); every shipped config has `χ=0.0` so this
is a documented no-op, not a silent numeric loss.
"""
function compiled_event_rates!(
    alpha, pi_, cols,
    S, E, I, R;
    β, σ, γ, ω, ψ, χ, pop, model::MGPModel,
    _...,
)
    ellE, ellI = ell(cols)
    @assert I ≥ ellI && E ≥ ellE

    infection   = model.events[findfirst(e -> e.name == :infection, model.events)]
    progression = model.events[findfirst(e -> e.name == :progression, model.events)]
    recovery    = model.events[findfirst(e -> e.name == :recovery, model.events)]
    waning      = model.events[findfirst(e -> e.name == :waning, model.events)]

    x = (S = S, E = E, I = I, R = R)
    θ = (β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, χ = χ, N = pop)

    α_inf = Float64(infection.hazard(x, θ))
    alpha[2] = alpha[1] = α_inf
    α_prog = Float64(progression.hazard(x, θ))
    alpha[4] = alpha[3] = α_prog
    alpha[5] = @indicator(I > ellI, γ * (I - ellI))
    alpha[6] = Float64(waning.hazard(x, θ))

    # NaiveProposal pi_u -- "source/parent-deme" split (mers_filter_suite.tex
    # lines 880-885's SEIRS transmission analogue: pi^emptyset = (n-ell)/n,
    # pi^b = 1/n per tracked lineage, summing to ell/n over the branch).
    pi_[1] = @indicator(I > 0, 1 - ellI / I)   # infection, untracked-I parent (identity)
    pi_[2] = @indicator(I > 0, ellI / I)       # infection, tracked-I parent (cross)
    pi_[3] = @indicator(E > 0, 1 - ellE / E)   # progression, untracked-E parent (identity)
    pi_[4] = @indicator(E > 0, ellE / E)       # progression, tracked-E parent (cross)
    pi_[6] = pi_[5] = 1.0

    ℓvec = [ellE, ellI]
    nvec = [E, I]
    compiled_decay(model, x, θ, ℓvec, nvec)
end

"""
    compiled_regular_part!(cols, ll, t, dt, S, E, I, R; model, kwargs...) -> (ll, S, E, I, R)

Drop-in structural replacement for `NaiveSEIR.regular_part!`
(`seir_naive.jl:123-186`): IDENTICAL control flow and RNG-call order
(`rcateg`, `rand()`, and — inside the `k==2`/`k==4` branches — `rand(cols[...])`)
so that, for a FIXED seed, this function and `NaiveSEIR.regular_part!` draw
the same random numbers and hence the same trajectory. Every likelihood term
is instead computed via M01-M07's IR:
  - `k∈{1,3}` (identity/no-op): `apply_move! = log(Φ_identity)` (`q=1`, no
    lineage draw needed — `full_transitions`+`reduce_event_indicator`'s
    `:noop` `ReducedTransition`, per this file's header derivation, verified
    to algebraically collapse to `NaiveSEIR`'s own `log(1-ellE/E)` /
    `log(1-ellI/I)` inner terms).
  - `k∈{2,4}` (cross-deme): `apply_move! = log(Φ_cross) - log(q)`, `q = 1/ℓ`
    (the uniform per-lineage draw `seir_naive.jl:152`/`:163` also performs),
    matching `test/kli_proposal_test.jl`'s already-verified
    `boost(Φ_cross,1/n_I)==naive_factor` cross-check (M07 Part B) when
    combined with the shared `-log(pi[k])` line below.
  - `k==5` (recovery): `ll -= log(1-ellI/I)` kept as-is — this is
    `NaiveProposal`'s own within-branch bookkeeping (no `y`-changing
    coloring op for a DEATH event; not a `Φ_u`/boost quantity — DEATH events
    are, per M03/M07, decay/rate-governed, not saturation-governed), so it is
    reproduced directly rather than routed through `full_transitions`.
  - `k==6` (waning): no coloring effect, matches `NaiveSEIR` exactly.
Decay is `compiled_decay` (this file), resolving M07's flagged discrepancy.
"""
function compiled_regular_part!(
    cols, ll,
    t, dt,
    S, E, I, R;
    model::MGPModel,
    kwargs...,
)
    tf = t + dt
    if t < tf
        alpha = similar(Vector{Prob}, 6)
        pi_ = similar(Vector{Prob}, 6)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellE, ellI = ell(cols)
        while t < tf
            decay = compiled_event_rates!(
                alpha, pi_, cols,
                S, E, I, R;
                model = model, kwargs...,
            )
            k, s = rcateg(alpha .* pi_)
            step = -log(rand()) / s
            if k > 0 && t + step < tf
                ll -= decay * step + log(pi_[k])
                infection = model.events[findfirst(e -> e.name == :infection, model.events)]
                progression = model.events[findfirst(e -> e.name == :progression, model.events)]
                if k == 1
                    S -= 1
                    E += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    # Use the UNREDUCED IdentityTransition directly, NOT
                    # reduce_event_indicator's collapsed `:noop` group — the
                    # latter also sums in InlineSameDemeTransition's mass
                    # (the "tracked I-parent stays tracked I" saturation),
                    # which M03 documented is NEVER reachable during a
                    # REGULAR step (`test/kli_full_transitions_test.jl`
                    # "SEIR infection" comment: "seir_naive.jl's
                    # regular_part! never visits the InlineSameDeme or Fork
                    # saturations during a REGULAR step" — a Q_u
                    # regular-vs-singular gating fact M03 explicitly scoped
                    # out and left for this milestone). Collapsing them here
                    # would silently double-count mass naive's proposal
                    # structurally never proposes.
                    ts = full_transitions(infection, ℓpost, npost)
                    Φid = only(filter(t -> t isa IdentityTransition, ts)).phi
                    # Φ_identity = (I-side factor, == pi_[1] exactly, since
                    # infection's r has a nontrivial slot in BOTH demes) *
                    # (E-side factor). The outer line above already charges
                    # -log(pi_[1]) once; dividing it back out here leaves
                    # exactly the E-side factor `1-ellE/E_post`, matching
                    # `seir_naive.jl:149`'s own inner term bit-for-bit
                    # (verified: `test/kli_full_transitions_test.jl`'s "SEIR
                    # infection" comment, lines 228-230).
                    ll += log(Float64(Φid)) - log(pi_[1])
                elseif k == 2
                    ellI_pre = ellI
                    b = rand(cols[NaiveSEIR.Infec])
                    ellE, ellI = swap!(cols, NaiveSEIR.Infec, NaiveSEIR.Expos, b)
                    S -= 1
                    E += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    # CrossDemeTransition is a singleton per event/state (M04:
                    # "never collapsed with the noop group"), so using it
                    # directly vs. via reduce_event_indicator's `:cross`
                    # group gives the identical value; kept direct for
                    # symmetry with the k==1 fix above.
                    ts = full_transitions(infection, ℓpost, npost)
                    Φcr = only(filter(t -> t isa CrossDemeTransition, ts)).phi
                    ll += log(Float64(Φcr)) - log(1 / ellI_pre)
                elseif k == 3
                    E -= 1
                    I += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    # progression has r=[0,1] (no slot in its own ancestral
                    # deme E), so it never produces an InlineSameDemeTransition
                    # at all -- IdentityTransition alone here (no collapsing
                    # issue), kept direct for the same reason as k==1.
                    ts = full_transitions(progression, ℓpost, npost)
                    Φid = only(filter(t -> t isa IdentityTransition, ts)).phi
                    ll += log(Float64(Φid))
                elseif k == 4
                    ellE_pre = ellE
                    b = rand(cols[NaiveSEIR.Expos])
                    ellE, ellI = swap!(cols, NaiveSEIR.Expos, NaiveSEIR.Infec, b)
                    E -= 1
                    I += 1
                    ℓpost = [ellE, ellI]
                    npost = [E, I]
                    ts = full_transitions(progression, ℓpost, npost)
                    Φcr = only(filter(t -> t isa CrossDemeTransition, ts)).phi
                    ll += log(Float64(Φcr)) - log(1 / ellE_pre)
                elseif k == 5
                    ll -= log(1 - ellI / I)
                    I -= 1
                    R += 1
                elseif k == 6
                    R -= 1
                    S += 1
                end
                t += step
            else
                step = tf - t
                ll -= decay * step
                break
            end
        end
        @assert I ≥ ellI && E ≥ ellE
    end
    ll, S, E, I, R
end

"""
    compiled_filter_pomp(gen; β, σ, γ, ω, ψ, χ, pop, S0, E0, I0, R0)

Structural mirror of `NaiveSEIR.filter_pomp` (`seir_naive.jl:195-253`): SAME
`rinit`, SAME `logdmeasure`, SAME singular update (`NaiveSEIR.singular_part!`,
reused verbatim — see this file's header "SCOPING DECISION"), but
`rprocess`'s regular step calls `compiled_regular_part!` (this file) instead
of `NaiveSEIR.regular_part!`. `model` (a `PhyloPOMP.SEIR` `MGPModel`) is
threaded through as an extra `args`/kwarg so `compiled_regular_part!` can
look up `Event`s by name.
"""
compiled_filter_pomp(
    gen::Genealogy;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02, χ = 0.0,
    pop = 100,
    S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ), χ = Float64(χ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; S0, E0, I0, R0, pop, _...)
            m = pop/(S0+E0+I0+R0)
            (
                node = one(Name),
                ll = zero(Prob),
                cols = Coloring(NaiveSEIR.Demes),
                S = round(Int64, m*Float64(S0)),
                E = round(Int64, m*Float64(E0)),
                I = round(Int64, m*Float64(I0)),
                R = round(Int64, m*Float64(R0)),
                live = true,
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, geneal,
                t, dt,
                S, E, I, R, live,
                args...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                ll, S, E, I, R, live = NaiveSEIR.singular_part!(
                    cols, geneal, node, ll, live,
                    S, E, I, R;
                    args...,
                )
                if live && dt > 0 && isfinite(ll)
                    ll, S, E, I, R = compiled_regular_part!(
                        cols, ll, t, dt,
                        S, E, I, R;
                        model = SEIR, args...,
                    )
                end
                (; node = node+1, ll = ll, cols = cols,
                 S = S, E = E, I = I, R = R, live = live)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end
