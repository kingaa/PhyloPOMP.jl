# Birth-death-exposed-infectious model, as R phylopomp's BDEI. Demes: E (exposed), I (infectious).
# An infectious host infects at rate λ; the new host is exposed. Sampling removes the host.

@mgp BDEI begin
    compartments = (E, I)
    demes = (E, I)
    params = (σ, λ, μ, χ)

    @event progression rate=σ*E pop=(E=-1, I=+1) move=swap(E => I)
    @event birth       rate=λ*I pop=(E=+1)       move=fork(I => E, I)
    @event death       rate=μ*I pop=(I=-1)       move=chop(I)
    @event sampling    rate=χ*I pop=(I=-1)       move=sample_remove(I) kind=singular
end
