# SIR model declaration (the simplest `@mgp` example).
#
# Compartments: S (susceptible), I (infectious), R (recovered).
# Demes: I (only infectious individuals carry lineages).
#
# Jump marks (3 total):
#   infection — BIRTH, I parent sires an I child; r=(2)
#   recovery  — DEATH of an I lineage
#   sampling  — SAMPLE of an I lineage; non-destructive (pop=()), so the
#               sampled individual stays infectious (serial sampling)
#
# There is no hand-coded filter for this model; it exists to show that the
# generic simulator (src/simulate.jl) needs nothing but this event table.

@mgp SIR begin
    compartments = (S, I, R)
    demes = (I,)
    params = (β, γ, ψ, N)

    @event infection rate=β*S*I/N pop=(S=-1, I=+1) move=fork(I => I, I) kind=regular
    @event recovery  rate=γ*I     pop=(I=-1, R=+1) move=chop(I)         kind=regular
    @event sampling  rate=ψ*I     pop=()           move=sample(I)       kind=singular
end
