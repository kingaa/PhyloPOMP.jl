using PhyloPOMP, Random
using PhyloPOMP: MGPModel, Event, BIRTH, MIGRATION, DEATH, SAMPLE, NEUTRAL, @mgp, @demes

# (b) MTBD declared with @mgp: no susceptibles, every rate per infected host.
# The full model, with r1 and r2, is now PhyloPOMP.MTBD (src/examples/mgp_mtbd.jl).
@mgp MTBDprobe begin
    compartments = (I_1, I_2)
    demes = (I_1, I_2)
    params = (λ11, λ12, λ21, λ22, μ1, μ2, ψ1, ψ2, m12, m21)
    @event b11 rate=λ11*I_1 pop=(I_1=+1) move=fork(I_1 => I_1, I_1)
    @event b12 rate=λ12*I_1 pop=(I_2=+1) move=fork(I_1 => I_1, I_2)
    @event b21 rate=λ21*I_2 pop=(I_1=+1) move=fork(I_2 => I_2, I_1)
    @event b22 rate=λ22*I_2 pop=(I_2=+1) move=fork(I_2 => I_2, I_2)
    @event d1  rate=μ1*I_1  pop=(I_1=-1) move=chop(I_1)
    @event d2  rate=μ2*I_2  pop=(I_2=-1) move=chop(I_2)
    @event s1  rate=ψ1*I_1  pop=(I_1=-1) move=sample_remove(I_1) kind=singular
    @event s2  rate=ψ2*I_2  pop=(I_2=-1) move=sample_remove(I_2) kind=singular
    @event m12 rate=m12*I_1 pop=(I_1=-1, I_2=+1) move=swap(I_1 => I_2)
    @event m21 rate=m21*I_2 pop=(I_2=-1, I_1=+1) move=swap(I_2 => I_1)
end
θ = (λ11=1.2, λ12=0.3, λ21=0.2, λ22=0.9, μ1=0.5, μ2=0.5, ψ1=0.3, ψ2=0.3, m12=0.0, m21=0.0)
G = simulate(MTBDprobe, θ; x0=(I_1=1, I_2=0), graft=[1,0], tmax=6.0, rng=MersenneTwister(3))
println("MTBD via @mgp: nsample = ", G.nsample, ", nodes = ", length(G))

# (a1) hand-built BIRTH with sum(r)=1 and consistent Δ: accepted silently?
@demes BadDemes I
bad1 = MGPModel(:Bad1, [:S,:I], [:I], [
    Event(:half_birth, [:S=>-1], (x,θ)->0.5*x.I*x.S/10, [1], BIRTH, 1, Int[], true, false),
    Event(:samp, Pair{Symbol,Int}[], (x,θ)->0.5*x.I, [1], SAMPLE, 1, Int[], false, true),
])
try
    G1 = simulate(bad1, (;); x0=(S=10, I=1), graft=[1], tmax=5.0, rng=MersenneTwister(1))
    println("sum(r)=1 BIRTH accepted: nsample=", G1.nsample, " nodes=", length(G1))
catch e
    println("sum(r)=1 BIRTH rejected: ", sprint(showerror, e)[1:min(end,150)])
end

# (a2) Δ key typo: silently dropped?
bad2 = MGPModel(:Bad2, [:S,:I,:R], [:I], [
    Event(:rec, [:I=>-1, :Rr=>+1], (x,θ)->1.0*x.I, [0], DEATH, 1, Int[], true, false),
])
try
    tr = simulate_trajectory(bad2, (;); x0=(S=0, I=3, R=0), graft=[3], tmax=50.0, rng=MersenneTwister(1))
    println("typo Δ key accepted: final state = ", tr.states[end])
catch e
    println("typo Δ key rejected: ", sprint(showerror, e)[1:min(end,150)])
end

# (a3) SAMPLE with r[from]=2: what happens?
bad3 = MGPModel(:Bad3, [:I], [:I], [
    Event(:s, Pair{Symbol,Int}[], (x,θ)->1.0*x.I, [2], SAMPLE, 1, Int[], false, true),
    Event(:d, [:I=>-1], (x,θ)->1.0*x.I, [0], DEATH, 1, Int[], true, false),
])
try
    G3 = simulate(bad3, (;); x0=(I=2,), graft=[2], tmax=5.0, rng=MersenneTwister(1))
    println("SAMPLE r[from]=2 accepted: nsample=", G3.nsample)
catch e
    println("SAMPLE r[from]=2 rejected: ", sprint(showerror, e)[1:min(end,150)])
end

# (a4) from out of range
bad4 = MGPModel(:Bad4, [:I], [:I], [
    Event(:d, [:I=>-1], (x,θ)->1.0*x.I, [0], DEATH, 2, Int[], true, false),
])
try
    simulate(bad4, (;); x0=(I=2,), graft=[2], tmax=5.0, rng=MersenneTwister(1))
    println("from=2 accepted")
catch e
    println("from out of range: ", typeof(e), " ", sprint(showerror, e)[1:min(end,120)])
end

# (e) hazard call cost: dynamic dispatch through Event.hazard::Function
println("isconcretetype(fieldtype(Event,:hazard)) = ", isconcretetype(fieldtype(Event,:hazard)))
