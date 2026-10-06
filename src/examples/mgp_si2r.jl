# SI2R superspreading model. Demes: I_L (low-rate spreader), I_H (super-spreader).
# A host is sampled at rate ψ and removed with probability r.
# r = 0 is the model of si2r_model.qmd; r = 1 is R phylopomp's runSI2R (chi = ψ).

@mgp SI2R begin
    compartments = (S, I_L, I_H, R)
    demes = (I_L, I_H)
    params = (β, κ, γ, ω, ψ, η_L, η_H, N, r)

    @event TL rate=β*S*I_L/N  pop=(S=-1, I_L=+1) move=fork(I_L => I_L, I_L) kind=regular
    @event TH rate=κ*β*S*I_H/N  pop=(S=-1, I_L=+1) move=fork(I_H => I_L, I_H) kind=regular
    @event L rate=η_L*I_L  pop=(I_L=-1, I_H=+1) move=swap(I_L => I_H) kind=regular
    @event H rate=η_H*I_H  pop=(I_H=-1, I_L=+1) move=swap(I_H => I_L) kind=regular
    @event RL rate=γ*I_L  pop=(I_L=-1, R=+1) move=chop(I_L) kind=regular
    @event RH rate=γ*I_H  pop=(I_H=-1, R=+1) move=chop(I_H) kind=regular
    @event W rate=ω*R  pop=(R=-1, S=+1) move=none kind=regular
    @event SL rate=(1-r)*ψ*I_L  pop=() move=sample(I_L) kind=singular
    @event SH rate=(1-r)*ψ*I_H  pop=() move=sample(I_H) kind=singular
    @event SL_remove rate=r*ψ*I_L  pop=(I_L=-1) move=sample_remove(I_L) kind=singular
    @event SH_remove rate=r*ψ*I_H  pop=(I_H=-1) move=sample_remove(I_H) kind=singular
end
