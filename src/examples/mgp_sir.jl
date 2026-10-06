# SIR: one deme (I). Sampling is non-destructive (pop=()), so the host stays infectious.
# No hand-coded filter; src/simulate.jl runs it from the event table alone.

@mgp SIR begin
    compartments = (S, I, R)
    demes = (I,)
    params = (β, γ, ψ, N)

    @event infection rate=β*S*I/N pop=(S=-1, I=+1) move=fork(I => I, I) kind=regular
    @event recovery  rate=γ*I     pop=(I=-1, R=+1) move=chop(I)         kind=regular
    @event sampling  rate=ψ*I     pop=()           move=sample(I)       kind=singular
end
