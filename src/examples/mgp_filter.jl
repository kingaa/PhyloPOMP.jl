# mgp_filter.jl
# =============================================================================
# Generic KLI filter over an MGPModel (see mgp.jl).
#
#   The LOOP STRUCTURE (regular_step!, the singular dispatch, the population
#   updates) is generic and model-independent. The LIKELIHOOD MATH lives in
#   exactly the handful of "slot" functions below. Each slot's docstring cites
#   the King–Lin–Ionides equation it must implement.
#
# STATUS: compiler-checked scaffold. The weight-bearing slots — `kli_select`,
#   `kli_decay`, `apply_move!`, `singular_update!` — are UNVERIFIED STUBS that
#   `error(...)`. They are deliberately not filled in: their bodies are the
#   thing that must be derived from the paper and checked, not guessed. Filling
#   a stub does NOT make it correct until it is (a) verified term-by-term
#   against the cited equation and (b) cross-validated (see "Validation gates").
#
# Depends on the generic coloring layer already in PhyloPOMP.jl:
#   Coloring, swap!, chop!, fork!, plant!, ell   (src/coloring.jl)
#   rcateg                                        (src/rcateg.jl)
# and treats the state `x` as a NamedTuple over compartments (see mgp.jl).
#
# Primary source: StructuredMGPs.pdf. All equation numbers are to that paper.
# =============================================================================

# -----------------------------------------------------------------------------
# Mechanical helpers (no KLI algebra; safe to implement directly).
# -----------------------------------------------------------------------------

"""
    kli_hazard(ev, x, θ) -> Float64

Population hazard αᵤ(t,x) for event `ev` (KLI §2.5). This is one factor of the
filter driver βᵤ = αᵤ·πᵤ (Theorem 5; SEIRS driver assembled in Eq. 45). Purely
mechanical: it just evaluates the event's stored hazard closure.
"""
kli_hazard(ev::Event, x, θ) = Float64(ev.hazard(x, θ))

"""
    apply_pop(x, ev) -> NamedTuple

Apply the population stoichiometry `ev.Δ` to the state. Mechanical; no genealogy
effect. Returns a new state (the NamedTuple state is immutable).
"""
function apply_pop(x::NamedTuple, ev::Event)
    isempty(ev.Δ) && return x
    names = keys(x)
    delta = Dict(ev.Δ)
    values = map(name -> getfield(x, name) + get(delta, name, 0), names)
    NamedTuple{names}(values)
end

# =============================================================================
# WEIGHT-BEARING SLOTS — the KLI reduction. UNVERIFIED STUBS.
# Fill and verify each against the cited equation before use.
# =============================================================================
"""
    kli_select(ev, cols, x) -> Float64

Return the KLI selection factor pi_u for regular event `ev`. Together with the
stored population hazard alpha_u, this gives the driver rate beta_u = alpha_u*pi_u
in Eq. 45. For SEIR this slot contains the lineage-count factors represented by
`pi` in `seir_naive.jl`; for example recovery gives `(I-ellI)/I`. Sampling is
singular and therefore has no regular driver entry.

The conditional choice among compatible coloring outcomes remains in
`apply_move!`. This separation is the one used by the technical companion.

UNVERIFIED STUB.
"""
function kli_select(ev::Event, cols, x)
    error("kli_select: unverified stub — implement per KLI Eqs. 44–46")
end

"""
    kli_decay(alpha, pi, cols, x, model) -> Float64

Total decay rate λ(t,x) accumulated between genealogy events (subtracted as
λ·Δt over each interval). Collects the exit hazards of events whose unobserved
occurrence would be inconsistent with the observed tree.

Reference: SEIRS decay KLI Eq. 47 — λ = ∫αSample dx' + ∫αRecov·1_{I≤ℓI} dx';
general definition in Appendix B, Eq. B2. Absorbs the return value of
`event_rates!` and the `ll -= decay*step` term in seir_naive.jl.

UNVERIFIED STUB.
"""
function kli_decay(alpha, pi, cols, x, model::MGPModel)
    error("kli_decay: unverified stub — implement per KLI Eq. 47 / Eq. B2")
end

"""
    apply_move!(cols, ev, x) -> Float64

Apply the coloring move implied by `ev` to the pruned genealogy `cols`. It
samples a compatible coloring outcome with conditional probability `q`, mutates
`cols` using `swap!`, `chop!`, or `fork!`, and returns
`log(phi_u) - log(q)`. The caller separately contributes `-log(pi_u)` for the
event-selection factor returned by `kli_select`.

For a naive proposal, `q` contains the uniform lineage or orientation choices
(`1/I`, `1/E`, or `1/2`). This is where the lineage multiplier on a regular
cross-deme swap and the factor two on a symmetric singular fork enter. The sum
of this return value and the caller term is the full importance correction.

Reference: boost B = Ψᵤ = ϕᵤ/πᵤ (Theorem 5); ϕᵤ = (binomial ratio)·Qᵤ (Eq. 9);
binomial ratio (Eq. 7) with the Chu–Vandermonde identity (Eq. 8); SEIRS
binomials (Eq. 37); compatibility Qᵤ (Eq. 38); assembled SEIRS boost (Eq. 46).
This is where the two recurring MERS defects belong:
  • the +log(ℓ) owed on a regular cross-deme swap (the ℓ-multiplier in the
    summed-over-m transmission term, Eq. 43, collapsed by Eq. 8);
  • the +log(2) owed on a singular symmetric fork (two-orientation branch-point
    boost, Eq. 41; equivalently the (ⁿ ℓ; 2 2)=2/(I(I-1)) ratio, Eq. 37).

UNVERIFIED STUB.
"""
function apply_move!(cols, ev::Event, x)
    error("apply_move!: unverified stub — implement per KLI Eqs. 9, 37, 38, 46")
end

"""
    singular_update!(cols, node, x, θ, model) -> (Δll, x′)

Singular update at an observed genealogy event (KLI Eq. 14; general Eq. B3).
Dispatch on the node type and apply the corresponding boost:
    branch-point → fork κ,   sample → chop χ,   inline → swap σ   (§5.1).
Roots are NOT handled here: they enter through the filter's initial condition
w(0,x,y) = p₀(x)·1{q(x,y)>0} (Theorem 5).

Reference: general singular part Eq. 14; worked SEIRS singular Eqs. 39–41;
Appendix B singular form Eq. B3. Absorbs the `singular_part!` boosts in
seir_naive.jl (log(β·S·I/pop) … -log(E·I); log(ψ·I)+log(1-ellI/I); log(χ·I)).

UNVERIFIED STUB.
"""
function singular_update!(cols, node, x, θ, model::MGPModel)
    error("singular_update!: unverified stub — implement per KLI Eqs. 14, 39–41, B3")
end

# =============================================================================
# GENERIC REGULAR STEP — structural loop (generic, safe). Calls the stubs above.
# Corresponds to KLI Eq. 13 (regular part of Theorem 5) / Eq. B2, integrated by
# the SMC scheme of Lemma 3 / Algorithm B1.
# =============================================================================

"""
    regular_step!(cols, ll, t, dt, x, model, θ) -> (ll, x, t)

Advance the filter across one inter-event interval [t, t+dt) with no observed
genealogy events. This is the forward SMC filter (Theorem 5, Eqs. 13–14,
integrated by Algorithm B1).

The naive / guided / hard proposals differ only in the importance kernel π
returned by `kli_select` and the conditional coloring proposal in `apply_move!`. The guided /
`relhaz` track builds an anticipatory π from a reverse-time sweep over the deme
process, to "borrow information from future events" (text p. 27; Discussion §6).
It is still this forward filter; only the proposal changes. It is NOT the adjoint.

(The adjoint filter — Corollary 6, Eqs. 15–16 / B11–B12 — is a separate BACKWARD
reformulation that computes the same L by integrating F from F⁻(T,x)=1. The
forward SMC scheme here does not use it. Reverse-time in `relhaz` refers to the
construction of π, not to the adjoint.)
"""
function regular_step!(cols, ll::Float64, t::Float64, dt::Float64, x, model::MGPModel, θ)
    tf = t + dt
    while t < tf
        alpha = Float64[ev.regular ? kli_hazard(ev, x, θ) : 0.0 for ev in model.events]
        pi = Float64[ev.regular ? kli_select(ev, cols, x) : 0.0 for ev in model.events]
        rates = alpha .* pi
        decay = kli_decay(alpha, pi, cols, x, model)
        total = sum(rates)
        if total <= 0
            ll -= decay * (tf - t)
            break
        end
        step = -log(rand()) / total
        if t + step < tf
            k, _ = rcateg(rates)
            ll -= decay*step + log(pi[k])
            x = apply_pop(x, model.events[k])
            ll += apply_move!(cols, model.events[k], x)
            t += step
        else
            ll -= decay * (tf - t)
            break
        end
    end
    ll, x, t
end

# =============================================================================
# Validation gates (from the KLI paper and the project's testing discipline).
# A filled-in filter is not trusted until it passes ALL of:
#   1. Known special cases reduce correctly:
#        • Kingman coalescent / Moran   — §4.3.1, Eqs. 17–20
#        • linear birth–death (Stadler) — §4.3.2, Eqs. 21–26
#   2. SIRS / SEIRS worked examples reproduced — SIRS Eqs. 32 & 36;
#      SEIRS singular Eqs. 39–41, regular Eqs. 42–43, driver/boost/decay 45–47.
#   3. Three-filter agreement (naive / guided / hard) on the same genealogy.
#   4. Chu–Vandermonde collapse (Eq. 8) verified where m/s sums are marginalized
#      — not left as permanently separate "inline-node inserted" states.
# =============================================================================
