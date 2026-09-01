# The two lineage coordinates are `(I_c, I_h)`. Their indices correspond to
# `(Camel, Human)` in the hand-written MERS filters.
@mgp MERS begin
    compartments = (S_c, I_c, S_h, I_h)
    demes = (I_c, I_h)
    params = (β_cc, β_ch, β_hc, β_hh, γ_c, γ_h, χ_c, χ_h,
              B_c, B_h, N_c, N_h)

    @event transmission_cc rate=β_cc*S_c*I_c/N_c pop=(S_c=-1, I_c=+1) move=fork(I_c => I_c, I_c) kind=regular
    @event transmission_hh rate=β_hh*S_h*I_h/N_h pop=(S_h=-1, I_h=+1) move=fork(I_h => I_h, I_h) kind=regular
    @event transmission_hc rate=β_hc*S_h*I_c/N_c pop=(S_h=-1, I_h=+1) move=fork(I_c => I_c, I_h) kind=regular
    @event transmission_ch rate=β_ch*S_c*I_h/N_h pop=(S_c=-1, I_c=+1) move=fork(I_h => I_c, I_h) kind=regular

    @event removal_c rate=γ_c*I_c pop=(I_c=-1) move=chop(I_c) kind=regular
    @event removal_h rate=γ_h*I_h pop=(I_h=-1) move=chop(I_h) kind=regular

    @event sampling_c rate=χ_c*I_c pop=(I_c=-1) move=sample_remove(I_c) kind=singular
    @event sampling_h rate=χ_h*I_h pop=(I_h=-1) move=sample_remove(I_h) kind=singular

    @event birth_c rate=B_c pop=(S_c=+1) move=none kind=regular
    @event birth_h rate=B_h pop=(S_h=+1) move=none kind=regular
    @event death_c rate=B_c*S_c/N_c pop=(S_c=-1) move=none kind=regular
    @event death_h rate=B_h*S_h/N_h pop=(S_h=-1) move=none kind=regular
end
