# mgp_mers_filter.jl
# =============================================================================
# M09: End-to-End MERS Compiler — applying M08's SEIR findings (decay-leftover,
# Q_u regular/singular gating) to MERS's two-deme (camel/human), 12-event
# table, and finding TWO NEW MERS-specific wrinkles along the way (documented
# below). Structural mirror of `mgp_seir_filter.jl` (M08): assembles M01-M07's
# pieces (`Event`/`MGPModel`, `full_transitions`/`reduce_event_indicator`
# (Phi_u), `total_decay`/`compiled_decay` (lambda), `driver`/`boost`) into an
# EXECUTABLE filter for MERS, built to be directly comparable against
# `mers_naive.jl` (read-only oracle, never modified).
#
#   SCOPING DECISION (same as M08, cited directly): this file compiles the
#   REGULAR part of the filter (`mers_compiled_regular_part!`, mirroring
#   `mers_naive.jl`'s `event_rates!`/`regular_part!`) via the compiler's own
#   machinery. The SINGULAR part (root-planting / sample-chop / branch-point
#   fork, `mers_naive.jl`'s `singular_part!`) is REUSED VERBATIM from
#   `NaiveMERS.singular_part!` rather than re-derived -- M08's justification
#   applies unchanged here: M02-M06 already established (Gates 2/3) that
#   `full_transitions`/`reduce_event_indicator` reproduce the same Phi_u
#   numbers `singular_part!` implicitly relies on; re-deriving the singular
#   branch-point/chop/root logic here would duplicate already-verified work
#   without new coverage. `mgp_filter.jl`'s generic stubs (`kli_select`,
#   `kli_decay`, `apply_move!`, `singular_update!`) remain untouched (still
#   `error()`) for the identical reason M08 gave: a canonical dispatch
#   signature for every EventType is out of scope for a per-model milestone.
#
# =============================================================================
# FINDING 1 (M08's decay-leftover pattern) -- CONFIRMED, NOT JUST ASSUMED,
# TO GENERALIZE DIRECTLY TO MERS'S TWO DEATH EVENTS
# =============================================================================
#
# `mgp_decay.jl`'s `total_decay`/`decay_contribution` and `mgp_seir_filter.jl`'s
# `compiled_decay` are BOTH already fully generic over `model.events` (neither
# is SEIR-specific despite living in the SEIR file) -- `compiled_decay` loops
# `for event in model.events; event.type == DEATH || continue; ...` with no
# SEIR-only assumption anywhere in its body. This milestone therefore does
# NOT duplicate `compiled_decay`'s code: `mers_compiled_event_rates!` below calls
# the SAME, already-exported `compiled_decay(model, x, θ, ℓ, n)` directly,
# passing MERS's own `x`/`θ`/`ℓ`/`n`. Verified (not assumed) by direct
# numerical reproduction of M07's own hand-derived MERS instance
# (`handoffs/M07_driver_boost_decay.md` Sec.1): gamma_c=1/2, gamma_h=2/5,
# chi_c=1/10, chi_h=3/10, I_C=5>ell_C=2, I_H=3=ell_H=3 (boundary):
#   total_decay(MERS,...)              = 13/5   (tex-literal lambda)
#   leftover(removal_c) = gamma_c*I_C*1{I_C>ell_C} - gamma_c*(I_C-ell_C)
#                        = 5/2 - gamma_c*3 = 5/2 - 3/2 = 1
#   leftover(removal_h) = 0 (I_H==ell_H, "above" gate false, both terms 0)
#   compiled_decay(MERS,...) = 13/5 + 1 + 0 = 18/5
# -- EXACTLY M07/M08's cited `mers_naive.jl` total (M08 handoff line 55: "13/5
# tex + 1 leftover(removal_c) + 0 leftover(removal_h) = 18/5"). Reproduced
# programmatically in `test/kli_mers_compiled_test.jl`'s "M07 decay instance"
# testset, in exact `Rational{Int}` arithmetic, not floating point.
#
# SAMPLE events (`sampling_c`/`sampling_h`) need NO leftover correction --
# confirmed against `mers_naive.jl:159` directly: `chi_c*Ic + chi_h*Ih` is
# added UNCONDITIONALLY, with no companion "reduced-rate" proposal anywhere in
# `event_rates!`/`regular_part!` for these marks (SAMPLE events are singular,
# `event.regular==false`, `mgp_filter.jl`'s regular driver already zeroes
# their `alpha`/`pi` -- they never enter the regular proposal at all, so there
# is no "regular proposal used a reduced rate" gap to leftover-compensate;
# `total_decay`'s `gate=1//1` SAMPLE branch already IS the naive `chi_d*I_d`
# term, term for term). This matches `mgp_decay.jl`'s own header (already
# vetted at M07) -- reconfirmed here specifically for MERS's `chi_c`/`chi_h`,
# not just asserted by analogy to SEIR's `psi*I`.
#
# =============================================================================
# FINDING 2 (M08's Q_u regular/singular-gating pattern) -- CONFIRMED FOR
# THC/TCH (structurally identical to SEIR's `infection`), BUT A GENUINELY NEW
# THIRD CASE FOUND FOR TCC/THH (the within-deme fork marks)
# =============================================================================
#
# `transmission_hc` (THC, r=(1,1), from=Camel) and `transmission_ch` (TCH,
# r=(1,1), from=Human) are STRUCTURALLY IDENTICAL to SEIR's `infection`
# (also r=(1,1), one slot in the ancestral/parent deme, one in the target
# deme) -- confirmed directly against `mers_naive.jl`'s actual code, not
# assumed from the shared `r` shape:
#   - `mers_naive.jl`'s k==3 (THC, untracked camel parent) / k==5 (TCH,
#     untracked human parent) branches perform NO lineage draw and use ONLY
#     `IdentityTransition` (s=(0,0)) -- `Phi_id/pi[k]` collapses to exactly
#     the naive inline term (verified: THC at ell_C=2,ell_H=1,I_C=5,I_H=4,
#     Phi_id=9/20, pi[3]=(1-2/5)=3/5, Phi_id/pi[3]=3/4 == naive's raw
#     `1-ell_H/I_H = 1-1/4 = 3/4`, exact match).
#   - `mers_naive.jl`'s k==4 (THC) / k==6 (TCH) branches perform EXACTLY ONE
#     lineage draw (`rand(cols[Camel])`/`rand(cols[Human])`) and realize
#     ONLY `CrossDemeTransition` (never `InlineSameDemeTransition`) --
#     `boost(Phi_cross, 1/ell_pre) == Phi_cross*ell_pre` matches the naive
#     `log(ell_pre) + log(Phi_cross)` combination exactly (verified: THC
#     instance above, `Phi_cross*ell_C = (3/20)*2 = 3/10`, matching the
#     `log(ell_C)+log(Phi_cross)` sum algebraically).
#   `InlineSameDemeTransition` (the "tracked parent, but the lineage
#   stays/continues in its OWN ancestral deme rather than transferring to the
#   new cross-deme child") is therefore NEVER realized by `mers_naive.jl`'s
#   regular proposal for THC/TCH -- exactly M08's SEIR finding, now confirmed
#   independently for MERS's cross-deme marks by direct inspection of
#   `mers_naive.jl`'s k==3/4/5/6 branches (not assumed to transfer from SEIR).
#
# `transmission_cc` (TCC, r=(2,0), BOTH slots in Camel) and `transmission_hh`
# (THH, r=(0,2), BOTH slots in Human) are DIFFERENT, and this is a genuinely
# NEW finding this milestone made (the task explicitly flagged this as
# unverified): `mers_naive.jl`'s k==1 (TCC) / k==2 (THH) branches perform NO
# lineage draw at all (`pi[1]=pi[2]=1`, a single undifferentiated branch) and
# their inline `ll +=` term is `1 - C(ell_d,2)/C(n_d,2)` -- this does NOT
# equal `Phi_noop` (M04's `reduce_event_indicator` collapsed
# `Identity+InlineSameDeme` group, a PLAIN SUM with no lineage-count
# weighting: verified numerically these DISAGREE, e.g. ell_C=2,I_C=5:
# `reduce_event_indicator` gives `Phi_noop=3/5`, but naive's raw term is
# `9/10`). The actual identity: because a TCC/THH event never distinguishes
# WHICH of the `ell_d` tracked lineages is silently "the same-deme parent"
# (the resulting `cols` is UNCHANGED regardless -- `IdentityTransition` and
# `InlineSameDemeTransition` are both `y'=y`), the aggregate probability of
# "no genealogy-visible fork" must be WEIGHTED by the number of ways each
# saturation's occupied slots could be filled by the `ell_d` tracked
# lineages -- `C(ell_d,0)=1` for Identity, `C(ell_d,1)=ell_d` for
# InlineSameDeme:
#     P_noop^aggregate = C(ell_d,0)*Phi_identity + C(ell_d,1)*Phi_inline
#                       = Phi_identity + ell_d*Phi_inline
# Verified to equal `mers_naive.jl`'s raw formula EXACTLY, at FOUR different
# `(ell_d,n_d)` instances (not just the one M07 instance), including a
# Chu-Vandermonde sum-to-1 cross-check (`Phi_id + ell_d*Phi_inline +
# C(ell_d,2)*Phi_fork == 1` exactly, confirming the general weighting
# convention, not a coincidence at one instance) -- see
# `test/kli_mers_compiled_test.jl`'s "TCC/THH weighted aggregate" testset.
# This is DISTINCT from THC/TCH's pattern (which needed only the standard
# `boost(Phi,1/ell)` -- itself already exactly `C(ell,1)*Phi`, so THC/TCH's
# cross case is actually the SAME `C(ell,s)*Phi_u(s)` weighting, just for a
# single split-out branch with its own explicit lineage draw rather than an
# aggregated, drawless sum): the general rule this milestone confirms is
# "naive's realized regular-time mass for a given target `y'` is
# `Sum_{s: y'(s)=y'} C(ell,s)*Phi_u(s)`, taken either as an aggregate sum (no
# lineage draw needed, when the group's outcome doesn't depend on WHICH
# lineages fill the slots -- TCC/THH's noop group) or as a single boosted
# term with an explicit uniform lineage draw (when it does -- THC/TCH's cross
# case)." `ForkTransition` (both slots tracked) is NEVER realized regularly
# by ANY of the four MERS BIRTH marks -- confirmed directly (no `mers_naive.jl`
# k-branch ever calls both `chop!`-adjacent bookkeeping AND a 2-lineage fork
# at regular time; `fork!` only appears in `singular_part!`, the SINGULAR/
# observed branch-point code, never in `regular_part!`).
#
# =============================================================================
# FINDING 3 (M09; RESOLVED post-M09, see below): `mgp_mers.jl` vs
# `mers_naive.jl` hazard-formula MISMATCH for the NEUTRAL `death_c`/`death_h`
# events -- originally reported and worked around, now fixed at the source
# =============================================================================
#
# As originally found at M09, `mgp_mers.jl:22-23` declared:
#     @event death_c rate=(S_c > 0 ? B_c : 0.0) pop=(S_c=-1) move=none
#     @event death_h rate=(S_h > 0 ? B_h : 0.0) pop=(S_h=-1) move=none
# i.e. a plain presence-gated CONSTANT rate `B_c`/`B_h`. But `mers_naive.jl`'s
# actual `event_rates!` (`mers_naive.jl:149-150`) computes
#     alpha[11] = Bc*Sc/Nc
#     alpha[12] = Bh*Sh/Nh
# a genuinely DIFFERENT (population-proportional) rate -- NOT the same
# formula (`S_c>0 ? B_c : 0.0` and `Bc*Sc/Nc` agree only in the degenerate
# case `Sc==Nc` for every `Sc>0`, never true in general). `birth_c`/`birth_h`
# (`mgp_mers.jl:20-21`, `rate=B_c`/`rate=B_h`) DID match `mers_naive.jl:147-148`
# (`alpha[9]=Bc`, `alpha[10]=Bh`) exactly -- only the two DEATH-side NEUTRAL
# marks disagreed. At M09 time this was a pre-existing gap in `mgp_mers.jl`
# (an earlier milestone's model-layer file, out of scope to modify at M09),
# so `mers_compiled_event_rates!`'s NEUTRAL-event alpha values (`alpha[11]`,
# `alpha[12]`) were computed directly from `mers_naive.jl`'s OWN formula
# (`Bc*Sc/Nc`, `Bh*Sh/Nh`) as a workaround, bypassing `model.events[...].hazard`
# for those two slots only.
#
# RESOLUTION: a later milestone fixed `mgp_mers.jl:22-23` itself to declare
# `rate=B_c*S_c/N_c`/`rate=B_h*S_h/N_h`, matching `mers_naive.jl`'s formula
# exactly (also documented as the "v28" fix in `mers_filter_suite.tex`'s
# changelog, which independently derived the same per-capita formula). With
# `MERS`'s own `death_c`/`death_h` hazards now correct, the workaround below
# is no longer needed: `mers_compiled_event_rates!` now sources `alpha[11]`/
# `alpha[12]` from `model.events[...].hazard` generically, exactly like every
# other event (including `birth_c`/`birth_h`). This was verified to preserve
# bit-exact Gate-5 numerical equivalence (`test/kli_mers_compiled_test.jl`) --
# `event.hazard` now produces the identical value the workaround used to
# hardcode, so the regular-step RNG stream is unaffected.
# =============================================================================

export mers_compiled_event_rates!, mers_compiled_regular_part!, mers_compiled_filter_pomp

"""
    _phi_of(ts::Vector{KLITransition}, ::Type{T}) -> Rational{Int}

Find the (at most one) `T <: KLITransition` in `ts` and return its `.phi`,
or `0//1` if `T` is absent from `ts` entirely. Needed because
`InlineSameDemeTransition`/`CrossDemeTransition`/`ForkTransition` are NOT
always enumerated by `full_transitions` -- `enumerate_saturations` collapses
a deme's feasible saturation range to `{0}` whenever `ell_d==0` for that
deme (`mgp_phi.jl`'s own documented boundary case), so e.g. TCC's
`InlineSameDemeTransition` (needs `s_C=1 <= ell_C`) legitimately does not
appear in `ts` when `ell_C==0` -- `only(filter(...))` would throw on an
empty collection in that case; the correct probability contribution is `0`,
not an error (`0` tracked lineages trivially cannot occupy 1 slot).
`IdentityTransition` (s all zeros) is always present regardless of `ell`, so
this helper is never actually needed for it, but is used uniformly below
for symmetry/robustness.
"""
function _phi_of(ts::AbstractVector{<:KLITransition}, ::Type{T}) where {T<:KLITransition}
    i = findfirst(t -> t isa T, ts)
    isnothing(i) ? zero(Rational{Int}) : ts[i].phi
end

"""
    mers_compiled_event_rates!(alpha, pi_, cols, Sc, Ic, Sh, Ih; Beta_cc, Beta_ch,
                           Beta_hc, Beta_hh, gamma_c, gamma_h, chi_c, chi_h,
                           Bc, Bh, Nc, Nh, model, _...) -> decay

Drop-in structural replacement for `NaiveMERS.event_rates!` (`mers_naive.jl:
132-162`), same 12-slot `alpha`/`pi` layout and semantics (1=TCC, 2=THH,
3/4=THC identity/cross, 5/6=TCH identity/cross, 7=removal_c, 8=removal_h,
9=birth_c, 10=birth_h, 11=death_c, 12=death_h). Every BIRTH/DEATH hazard is
sourced from `model.events[...].hazard` (`MERS`'s own `Event` table), including
`death_c`/`death_h` -- see "FINDING 3" above for the M09-era workaround this
used to require and how it was resolved by fixing `mgp_mers.jl` at the source.
The return value is `compiled_decay` (`mgp_seir_filter.jl`, already fully
generic -- see "FINDING 1" above for why no MERS-specific decay function is
defined in this file).
"""
function mers_compiled_event_rates!(
    alpha, pi_, cols,
    Sc, Ic, Sh, Ih;
    Beta_cc, Beta_ch, Beta_hc, Beta_hh,
    gamma_c, gamma_h, chi_c, chi_h, Bc, Bh, Nc, Nh,
    model::MGPModel,
    _...,
)
    ellc, ellh = ell(cols)
    @assert Ic ≥ ellc && Ih ≥ ellh

    tcc = model.events[findfirst(e -> e.name == :transmission_cc, model.events)]
    thh = model.events[findfirst(e -> e.name == :transmission_hh, model.events)]
    thc = model.events[findfirst(e -> e.name == :transmission_hc, model.events)]
    tch = model.events[findfirst(e -> e.name == :transmission_ch, model.events)]
    birth_c = model.events[findfirst(e -> e.name == :birth_c, model.events)]
    birth_h = model.events[findfirst(e -> e.name == :birth_h, model.events)]
    death_c = model.events[findfirst(e -> e.name == :death_c, model.events)]
    death_h = model.events[findfirst(e -> e.name == :death_h, model.events)]

    x = (S_c = Sc, I_c = Ic, S_h = Sh, I_h = Ih)
    θ = (β_cc = Beta_cc, β_ch = Beta_ch, β_hc = Beta_hc, β_hh = Beta_hh,
         γ_c = gamma_c, γ_h = gamma_h, χ_c = chi_c, χ_h = chi_h,
         B_c = Bc, B_h = Bh, N_c = Nc, N_h = Nh)

    alpha[1] = Float64(tcc.hazard(x, θ))
    alpha[2] = Float64(thh.hazard(x, θ))
    alpha[4] = alpha[3] = Float64(thc.hazard(x, θ))
    alpha[6] = alpha[5] = Float64(tch.hazard(x, θ))
    alpha[7] = @indicator(Ic > ellc, gamma_c * (Ic - ellc))
    alpha[8] = @indicator(Ih > ellh, gamma_h * (Ih - ellh))
    # FINDING 3 (resolved): mgp_mers.jl's death_c/death_h hazards were fixed
    # to match mers_naive.jl's formula (Bc*Sc/Nc, Bh*Sh/Nh), so all four
    # NEUTRAL events now go through model.events[...].hazard generically,
    # same as every other event.
    alpha[9]  = Float64(birth_c.hazard(x, θ))
    alpha[10] = Float64(birth_h.hazard(x, θ))
    alpha[11] = Float64(death_c.hazard(x, θ))
    alpha[12] = Float64(death_h.hazard(x, θ))

    pi_[1] = one(Prob)
    pi_[2] = one(Prob)
    pi_[3] = @indicator(Ic > 0, 1 - ellc / Ic)
    pi_[4] = @indicator(Ic > 0, ellc / Ic)
    pi_[5] = @indicator(Ih > 0, 1 - ellh / Ih)
    pi_[6] = @indicator(Ih > 0, ellh / Ih)
    for k in 7:12
        pi_[k] = one(Prob)
    end

    ℓvec = [ellc, ellh]
    nvec = [Ic, Ih]
    compiled_decay(model, x, θ, ℓvec, nvec)
end

"""
    mers_compiled_regular_part!(cols, ll, t, dt, Sc, Ic, Sh, Ih; model, kwargs...)
        -> (ll, Sc, Ic, Sh, Ih)

Drop-in structural replacement for `NaiveMERS.regular_part!` (`mers_naive.jl:
164-244`): IDENTICAL control flow (including the `while t+step < tf`/inner
`if t+step < tf`/`else` shape, NOT simplified to `while t < tf` --
`mers_naive.jl`'s own loop condition is mirrored exactly, since deviating
from it would change how many `rand()`-consuming draws occur and desync the
seed-matched RNG stream) and RNG-call order, so that, for a FIXED seed, this
function and `NaiveMERS.regular_part!` draw the same random numbers and hence
the same trajectory. Every likelihood term is instead computed via M01-M07's
IR:
  - `k∈{1,2}` (TCC/THH, within-deme fork marks): the `C(ell_d,s)`-weighted
    aggregate `Phi_identity + ell_d*Phi_inline` derived in this file's header
    "FINDING 2" -- a genuinely MERS-specific pattern, NOT the SEIR
    `boost(Phi,pi)` shape (there is no lineage draw to divide out here).
  - `k∈{3,5}` (THC/TCH identity): `log(Phi_identity) - log(pi_[k])`, the SAME
    divide-out pattern M08 used for SEIR's infection identity branch.
  - `k∈{4,6}` (THC/TCH cross): `log(Phi_cross) - log(1/ell_pre)`, the SAME
    `boost(Phi_cross, 1/ell_pre)` pattern M08 used for SEIR's infection cross
    branch, now confirmed for MERS's THC/TCH marks too.
  - `k∈{7,8}` (removal_c/removal_h, DEATH): kept as direct reproduction of
    naive's own bookkeeping (`ll -= log(1-ell_d/n_d)`) -- DEATH events are
    decay/rate-governed, not saturation-governed (M03/M07), so there is no
    `Phi_u`/boost term to route through `full_transitions` here, exactly as
    M08 did for SEIR's recovery.
  - `k∈{9,...,12}` (birth_c/birth_h/death_c/death_h, NEUTRAL): plain
    population updates, no coloring effect, no `ll` term at all -- matches
    `mers_naive.jl:223-230` exactly.
Decay is `compiled_decay` (already-generic, from `mgp_seir_filter.jl` -- see
"FINDING 1" above).
"""
function mers_compiled_regular_part!(
    cols, ll,
    t, dt,
    Sc, Ic, Sh, Ih;
    model::MGPModel,
    kwargs...,
)
    tf = t + dt
    if t < tf
        alpha = similar(Vector{Prob}, 12)
        pi_ = similar(Vector{Prob}, 12)
        step::Time = zero(Time)
        decay::Prob = zero(Prob)
        ellc, ellh = ell(cols)

        tcc = model.events[findfirst(e -> e.name == :transmission_cc, model.events)]
        thh = model.events[findfirst(e -> e.name == :transmission_hh, model.events)]
        thc = model.events[findfirst(e -> e.name == :transmission_hc, model.events)]
        tch = model.events[findfirst(e -> e.name == :transmission_ch, model.events)]

        while t + step < tf
            decay = mers_compiled_event_rates!(
                alpha, pi_, cols,
                Sc, Ic, Sh, Ih;
                model = model, kwargs...,
            )
            k, s = rcateg(alpha .* pi_)
            step = -log(rand()) / s
            if t + step < tf
                ll -= decay * step + log(pi_[k])
                if k == 1
                    Sc -= 1
                    Ic += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(tcc, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    Φinl = _phi_of(ts, InlineSameDemeTransition)
                    # C(ellc,0)*Φid + C(ellc,1)*Φinl = Φid + ellc*Φinl
                    # (FINDING 2: TCC's no-lineage-draw aggregate).
                    ll += log(Float64(Φid) + ellc * Float64(Φinl))
                elseif k == 2
                    Sh -= 1
                    Ih += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(thh, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    Φinl = _phi_of(ts, InlineSameDemeTransition)
                    ll += log(Float64(Φid) + ellh * Float64(Φinl))
                elseif k == 3
                    Sh -= 1
                    Ih += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(thc, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    ll += log(Float64(Φid)) - log(pi_[3])
                elseif k == 4
                    ellc_pre = ellc
                    b = rand(cols[NaiveMERS.Camel])
                    ellc, ellh = swap!(cols, NaiveMERS.Camel, NaiveMERS.Human, b)
                    Sh -= 1
                    Ih += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(thc, ℓpost, npost)
                    Φcr = _phi_of(ts, CrossDemeTransition)
                    ll += log(Float64(Φcr)) - log(1 / ellc_pre)
                elseif k == 5
                    Sc -= 1
                    Ic += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(tch, ℓpost, npost)
                    Φid = _phi_of(ts, IdentityTransition)
                    ll += log(Float64(Φid)) - log(pi_[5])
                elseif k == 6
                    ellh_pre = ellh
                    b = rand(cols[NaiveMERS.Human])
                    ellc, ellh = swap!(cols, NaiveMERS.Human, NaiveMERS.Camel, b)
                    Sc -= 1
                    Ic += 1
                    ℓpost = [ellc, ellh]
                    npost = [Ic, Ih]
                    ts = full_transitions(tch, ℓpost, npost)
                    Φcr = _phi_of(ts, CrossDemeTransition)
                    ll += log(Float64(Φcr)) - log(1 / ellh_pre)
                elseif k == 7
                    ll -= log(1 - ellc / Ic)
                    Ic -= 1
                elseif k == 8
                    ll -= log(1 - ellh / Ih)
                    Ih -= 1
                elseif k == 9
                    Sc += 1
                elseif k == 10
                    Sh += 1
                elseif k == 11
                    Sc -= 1
                elseif k == 12
                    Sh -= 1
                else
                    @assert false "impossible event" # COV_EXCL_LINE
                end
                t += step
            else
                step = tf - t
                ll -= decay * step
                break
            end
        end
        @assert Ic ≥ ellc && Ih ≥ ellh
    end
    ll, Sc, Ic, Sh, Ih
end

"""
    mers_compiled_filter_pomp(; Beta_cc, Beta_ch, Beta_hc, Beta_hh, gamma_c,
                          gamma_h, chi_c, chi_h, Bc, Bh, Sc0, Sh0, Ic0, Ih0,
                          Nc, Nh)

Structural mirror of `NaiveMERS.filter_pomp` (`mers_naive.jl:252-324`): SAME
`rinit`, SAME `logdmeasure`, SAME singular update (`NaiveMERS.singular_part!`,
reused verbatim -- see this file's header "SCOPING DECISION"), but
`rprocess`'s regular step calls `mers_compiled_regular_part!` (this file) instead
of `NaiveMERS.regular_part!`. `model` (`PhyloPOMP.MERS`, an `MGPModel`) is
threaded through as an extra kwarg so `mers_compiled_regular_part!` can look up
`Event`s by name. `gen` (the observed genealogy) is a required positional
argument, matching `mgp_seir_filter.jl`'s `mers_compiled_filter_pomp` signature
(unlike `NaiveMERS.filter_pomp`, which hardcodes `mers_tree` -- this keeps
`mers_compiled_filter_pomp` usable with simulated genealogies for Gate 5).
"""
mers_compiled_filter_pomp(
    gen::Genealogy;
    Beta_cc = 4.0, Beta_ch = 0.0, Beta_hc = 1.0, Beta_hh = 4.0,
    gamma_c = 1.0, gamma_h = 1.0,
    chi_c = 1.0, chi_h = 0.0,
    Bc = 0.1, Bh = 0.03,
    Sc0 = 1.0, Sh0 = 1.0,
    Ic0 = 0.01, Ih0 = 0.0,
    Nc = 10000, Nh = 10000,
) = begin
    pomp(
        params = (
            Beta_cc = Float64(Beta_cc), Beta_ch = Float64(Beta_ch),
            Beta_hc = Float64(Beta_hc), Beta_hh = Float64(Beta_hh),
            gamma_c = Float64(gamma_c), gamma_h = Float64(gamma_h),
            chi_c = Float64(chi_c), chi_h = Float64(chi_h),
            Bc = Float64(Bc), Bh = Float64(Bh),
            Sc0 = Float64(Sc0), Sh0 = Float64(Sh0),
            Ic0 = Float64(Ic0), Ih0 = Float64(Ih0),
            Nc = Float64(Nc), Nh = Float64(Nh),
        ),
        t0 = timezero(gen),
        times = times(gen),
        rinit = function (; Sc0, Sh0, Ic0, Ih0, Nc, Nh, _...)
            fc = Nc / (Sc0 + Ic0)
            fh = Nh / (Sh0 + Ih0)
            (
                node = one(Name),
                ll = zero(Prob),
                cols = Coloring(NaiveMERS.Demes),
                Sc = round(Int64, fc * Sc0),
                Ic = round(Int64, fc * Ic0),
                Sh = round(Int64, fh * Sh0),
                Ih = round(Int64, fh * Ih0),
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, geneal,
                Sc, Ic, Sh, Ih,
                t, dt,
                args...,
                )
                cols = copy(cols)
                ll = zero(Prob)
                ll, Sc, Ic, Sh, Ih = NaiveMERS.singular_part!(
                    cols, ll, geneal, node,
                    Sc, Ic, Sh, Ih;
                    args...,
                )
                if isfinite(ll)
                    ll, Sc, Ic, Sh, Ih = mers_compiled_regular_part!(
                        cols, ll, t, dt,
                        Sc, Ic, Sh, Ih;
                        model = MERS, args...,
                    )
                end
                (; node = node + one(Name), ll, cols, Sc, Ic, Sh, Ih)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (geneal = gen,),
    )
end
