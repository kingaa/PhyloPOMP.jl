# Two-type birth-death-sampling model (Kuehnert et al. 2016), as R phylopomp's MTBD2.
# No susceptibles: every rate is per infected host.
# A type-i host is sampled at rate psi_i and removed with probability r_i.
# Parameter names follow `epi_params` (mtbd_funs.jl), plus r1, r2.
# NaiveMTBD and GuidedMTBD handle only r1 = r2 = 1.

@mgp MTBD begin
    compartments = (I1, I2)
    demes = (I1, I2)
    params = (lambda11, lambda12, lambda21, lambda22, m12, m21, mu1, mu2, psi1, psi2, r1, r2)

    @event birth11 rate=lambda11*I1 pop=(I1=+1) move=fork(I1 => I1, I1)
    @event birth12 rate=lambda12*I1 pop=(I2=+1) move=fork(I1 => I1, I2)
    @event birth21 rate=lambda21*I2 pop=(I1=+1) move=fork(I2 => I2, I1)
    @event birth22 rate=lambda22*I2 pop=(I2=+1) move=fork(I2 => I2, I2)

    @event migrate12 rate=m12*I1 pop=(I1=-1, I2=+1) move=swap(I1 => I2)
    @event migrate21 rate=m21*I2 pop=(I2=-1, I1=+1) move=swap(I2 => I1)

    @event death1 rate=mu1*I1 pop=(I1=-1) move=chop(I1)
    @event death2 rate=mu2*I2 pop=(I2=-1) move=chop(I2)

    @event sample_remove1 rate=r1*psi1*I1     pop=(I1=-1) move=sample_remove(I1) kind=singular
    @event sample_remove2 rate=r2*psi2*I2     pop=(I2=-1) move=sample_remove(I2) kind=singular
    @event sample1        rate=(1-r1)*psi1*I1 pop=()      move=sample(I1)        kind=singular
    @event sample2        rate=(1-r2)*psi2*I2 pop=()      move=sample(I2)        kind=singular
end
