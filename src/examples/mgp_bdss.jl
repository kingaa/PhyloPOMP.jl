# Birth-death superspreading model, as R phylopomp's BDSS. Demes: N (normal), S (superspreader).
# Births keep the parent and add one host of the child's type. Sampling removes the host.

@mgp BDSS begin
    compartments = (N, S)
    demes = (N, S)
    params = (λ_nn, λ_ns, λ_sn, λ_ss, μ, χ)

    @event birth_nn   rate=λ_nn*N pop=(N=+1) move=fork(N => N, N)
    @event birth_ns   rate=λ_ns*N pop=(S=+1) move=fork(N => N, S)
    @event birth_sn   rate=λ_sn*S pop=(N=+1) move=fork(S => S, N)
    @event birth_ss   rate=λ_ss*S pop=(S=+1) move=fork(S => S, S)
    @event death_n    rate=μ*N    pop=(N=-1) move=chop(N)
    @event death_s    rate=μ*S    pop=(S=-1) move=chop(S)
    @event sampling_n rate=χ*N    pop=(N=-1) move=sample_remove(N) kind=singular
    @event sampling_s rate=χ*S    pop=(S=-1) move=sample_remove(S) kind=singular
end
