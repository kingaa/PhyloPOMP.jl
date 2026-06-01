module GuidedSEIRTest

import ..Main: h1, h2

@info h1("SEIR model with guided proposals")

using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using Test
using BenchmarkTools

@demes SEIR Expos Infec
using .SEIR: Expos, Infec, T as DemeType

using PhyloPOMP: Root, Node, Sample, Name, Time, repair!

seir_convert(n::GenealNode) = begin
    if n.type==Root
        n1 = GenealNode{DemeType}(n.name,n.slate,missing,n.type)
    elseif n.type==Node
        n1 = GenealNode{DemeType}(n.name,n.slate,Infec,n.type,n.parent)
    elseif n.type==Sample
        n1 = GenealNode{DemeType}(n.name,n.slate,Infec,n.type,n.parent)
    end
    n1.children = copy(n.children)
    n1
end

seir_convert(g::Genealogy) = begin
    h = Genealogy{DemeType}(g.t0,g.time)
    for n in g.nodes
        push!(h.nodes,seir_convert(n))
    end
    repair!(h)
    h
end

seir_singular(
    g, ll, cols;
    S, E, I, R, pop, β, ψ,
    _...,
) = begin
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    if g.type==Root
        if length(g.chillins) == 1
            if E-ellE+I-ellI > 0
                k, _, p = rcateg([E-ellE, I-ellI], true)
                ll -= log(p)
                if k==1
                    push!(cols[Expos], g.chillins[1])
                else
                    push!(cols[Infec], g.chillins[1])
                end
            else
                ll += Float64(-Inf)
                push!(cols[Infec], g.chillins[1])
                I += 1
            end
        else
            error("incompatible $(length(g.chillins)) > 1 at root $(g.name), t=$(g.time)")
        end
    elseif g.type==Sample
        if g.parlin ∉ cols[Infec]
            ll += Float64(-Inf)
            delete!(cols[Expos], g.parlin)
            push!(cols[Infec], g.parlin)
            E -= 1
            I += 1
        end
        ll += log(ψ*I)
        delete!(cols[Infec], g.parlin)
        if length(g.chillins) == 0
            ll += log(1-ell(cols,Infec)/I)
        elseif length(g.chillins) == 1
            push!(cols[Infec], g.chillins[1])
            ll -= log(I)
        else
            error("incompatible $(length(g.chillins)) > 1 at sample $(g.name), t=$(g.time)")
        end
    elseif g.type==Node
        if g.parlin ∉ cols[Infec]
            ll += Float64(-Inf)
            delete!(cols[Expos], g.parlin)
            push!(cols[Infec], g.parlin)
            E -= 1
            I += 1
        end
        if length(g.chillins) == 2
            ll += log(β*S*I/pop)
            delete!(cols[Infec], g.parlin)
            k, _, p = rcateg([g.present[1,1]*g.present[2,2], g.present[1,2]*g.present[2,1]], true)
            ll -= log(p)
            if k==1
                push!(cols[Expos], g.chillins[1])
                push!(cols[Infec], g.chillins[2])
            else
                push!(cols[Infec], g.chillins[1])
                push!(cols[Expos], g.chillins[2])
            end
            S -= 1
            E += 1
            ll -= log(E*I)
        else
            error("incompatible $(length(g.chillins)) ≠ 2 at node $(g.name), t=$(g.time)")
        end
    else
        @assert false "impossible node type"
    end
    (; ll = ll, cols = cols, S = S, E = E, I = I, R = R)
end

choice1(g, x, cols, rh) = begin
    lins = map(i->g.linmap[i],collect(cols))
    v = [x-ell(cols), sum(rh[lins])]
    v./sum(v)
end

choice2(g, cols, rh) = begin
    lins = map(i->g.linmap[i],collect(cols))
    k,_,p = rcateg(rh[lins],true)
    g.alllins[lins[k]], p
end

event_rates!(
    alpha, pi, cols, node, rh;
    S, E, I, R,
    β, σ, γ, ω, ψ, pop,
    _...,
) = begin
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    alpha[2] = alpha[1] = β*S*I/pop
    alpha[4] = alpha[3] = σ*E
    alpha[5] = @indicator(I > ellI, γ*I)
    alpha[6] = ω*R
    cI = choice1(node,I,cols[Infec],rh[(Infec,Expos)])
    pi[1] = @indicator(I > 0, cI[1])
    pi[2] = @indicator(I > 0, cI[2])
    cE = choice1(node,E,cols[Expos],rh[(Expos,Infec)])
    pi[3] = @indicator(E > 0, cE[1])
    pi[4] = @indicator(E > 0, cE[2])
    pi[6] = pi[5] = 1.0
    ψ*I + @indicator(I ≤ ellI, γ*I)
end

seir_regular(
    guide, node, ll, cols;
    S, E, I, R,
    β, σ, γ, ω, ψ, pop,
    _...,
) = begin
    alpha = similar(Vector{Float64}, 6)
    pi = similar(Vector{Float64}, 6)
    (ellE,ellI) = ell(cols)
    @assert I ≥ ellI && E ≥ ellE
    n = guide[node]
    t = n.tbeg
    tf = n.tend
    rh = relhaz(guide,t,node)
    if t < tf
        penalty = event_rates!(
            alpha, pi, cols, n, rh;
            S = S, E = E, I = I, R = R,
            β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
        )
        k, s = rcateg(alpha .* pi)
        step::Time = -log(rand())/s
        while (t+step < tf)
            ll -= penalty*step+log(pi[k])
            if k==1
                S -= 1
                E += 1
                ll += log(1-ellE/E)
            elseif k==2
                b,p = choice2(n,cols[Infec],rh[(Infec,Expos)])
                ll -= log(p)
                @assert b ∈ cols[Infec]
                delete!(cols[Infec], b)
                push!(cols[Expos], b)
                (ellE,ellI) = ell(cols)
                S -= 1
                E += 1
                ll += log(1-ellI/I)-log(E)
            elseif k==3
                E -= 1
                I += 1
                ll += log(1-ellI/I)
            elseif k==4
                b,p = choice2(n,cols[Expos],rh[(Expos,Infec)])
                ll -= log(p)
                @assert b ∈ cols[Expos]
                delete!(cols[Expos], b)
                push!(cols[Infec], b)
                (ellE,ellI) = ell(cols)
                E -= 1
                I += 1
                ll -= log(I)
            elseif k==5
                I -= 1
                R += 1
            elseif k==6
                R -= 1
                S += 1
            else
                @assert false "impossible error!"
            end
            @assert I ≥ ellI && E ≥ ellE
            t += step
            penalty = event_rates!(
                alpha, pi, cols, n, rh;
                S = S, E = E, I = I, R = R,
                β = β, σ = σ, γ = γ, ω = ω, ψ = ψ, pop = pop,
            )
            k, s = rcateg(alpha .* pi)
            step = -log(rand())/s
        end
        step = tf - t
        ll -= penalty*step
    end
    @assert I ≥ ellI && E ≥ ellE
    (; ll = ll, cols = cols, S = S, E = E, I = I, R = R)
end

seir(
    guide::Guide;
    β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02,
    pop = 100,
    S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
) = begin
    pomp(
        params = (
            β = Float64(β), σ = Float64(σ), γ = Float64(γ),
            ω = Float64(ω), ψ = Float64(ψ),
            pop = Float64(pop),
            S0 = Float64(S0), E0 = Float64(E0),
            I0 = Float64(I0), R0 = Float64(R0),
        ),
        t0 = timezero(guide),
        times = times(guide),
        rinit = function (; S0, E0, I0, R0, pop, _...)
            m = pop/(S0+E0+I0+R0)
            (
                node = one(Name),
                ll = zero(Float64),
                cols = Coloring(SEIR),
                S = round(Int64, m*Float64(S0)),
                E = round(Int64, m*Float64(E0)),
                I = round(Int64, m*Float64(I0)),
                R = round(Int64, m*Float64(R0)),
            )
        end,
        rprocess = onestep(
            function (
                ; node, ll, cols, guide,
                S, E, I, R,
                args...,
                )
                cols = copy(cols)
                ll = zero(Float64)
                ll, cols, S, E, I, R = seir_singular(
                    guide[node], ll, cols;
                    S = S, E = E, I = I, R = R,
                    args...,
                )
                ll, cols, S, E, I, R = seir_regular(
                    guide, node, ll, cols;
                    S = S, E = E, I = I, R = R,
                    args...,
                )
                (; node = node+1, ll = ll, cols = cols,
                 S = S, E = E, I = I, R = R)
            end,
        ),
        logdmeasure = function (; ll, _...)
            ll
        end,
        userdata = (guide = guide,),
    )
end

@testset verbose=true "SEIR model with guided proposals" begin

    g1 = parse_newick(readlines("seir1.nwk"), time = 50.0)
    @test g1 isa Genealogy{PhyloPOMP.Unstructured.T}

    g2 = seir_convert(g1)
    @test g2 isa Genealogy{DemeType}

    g3 = guide(g2,fsmarkov(Expos=>0.1,Infec=>1,(Expos,Infec)=>1))

    p = seir(g3)
    @test p isa POMP.PompObject

    @info h2("simulate test")
    sm = simulate(p, nsim = 3)
    @time sm = simulate(p, nsim = 3)
    @test sm isa Matrix{<:POMP.PompObject}

    @info h2("pfilter test")
    pf = pfilter(p, Np = 100)
    @time pf = pfilter(p, Np = 100)
    @test pf isa POMP.PfilterdPompObject
    @test isfinite(pf.logLik)

    @info h2("pfilter benchmark")
    @btime pfilter($p, Np = 1000)

    logmeanexp(x) = begin
        xmax = maximum(x)
        xmax + log(sum(exp.(x .- xmax))) - log(length(x))
    end

    using Statistics: std
    x = [logLik(pfilter(p,Np=1000)) for _ ∈ 1:10]
    @info "logLik = $(round(logmeanexp(x),digits=2)) ± $(round(std(x),sigdigits=3))"

end

end
