using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using Test

@info h1("SEIR model tests")

@testset verbose=true "SEIR model" begin

    g = parse_newick(readlines("seir1.nwk"),time=50.0)
    @test g isa Genealogy{PhyloPOMP.Unstructured.T}

    seir_singular(
        geneal,node;
        S,E,I,pop,β,ψ,linE,linI,ll,_...,
    ) = begin
        n = geneal[node]
        if n.type==PhyloPOMP.Root
            if E-length(linE)+I-length(linI) > 0
                k,_,p = rcateg([E-length(linE),I-length(linI)],true)
                ll -= log(p)
                if k==1
                    push!(linE,n.lineage)
                else
                    push!(linI,n.lineage)
                end
            else
                ll += Float64(-Inf)
                push!(linI,n.lineage)
                I += 1
            end
        elseif n.type==PhyloPOMP.Sample
            if n.lineage ∉ linI
                ll += Float64(-Inf)
                delete!(linE,n.lineage)
                push!(linI,n.lineage)
                E -= 1; I += 1;
            end
            ll += log(ψ*I)
            delete!(linI,n.lineage)
            @assert length(n.children)<2 "impossible: $(length(n.children)) > 1 at sample $(n.name), t=$(n.time)"
            if length(n.children) == 1
                push!(linI,geneal[n.children[1]].lineage)
                ll -= log(I)
            else
                ll += log(1-length(linI)/I)
            end
        elseif n.type==PhyloPOMP.Node
            if n.lineage ∉ linI
                ll += Float64(-Inf)
                delete!(linE,n.lineage)
                push!(linI,n.lineage)
                E -= 1; I += 1;
            end
            @assert length(n.children)==2 "impossible: $(length(n.children)) ≠ 2 at node $(n.name), t=$(n.time)"
            ll += log(β*S*I/pop)
            delete!(linI,n.lineage)
            k,_,p = rcateg([1,1],true)
            ll -= log(p)
            if k==1
                push!(linE,geneal[n.children[1]].lineage)
                push!(linI,geneal[n.children[2]].lineage)
            else
                push!(linI,geneal[n.children[1]].lineage)
                push!(linE,geneal[n.children[2]].lineage)
            end
            S -= 1; E += 1;
            ll -= log(E*I)
        end
        @assert (I ≥ length(linI) && E ≥ length(linE)) "one: ($I,$(length(linI)),$E,$(length(linI)))"
        (S=S,E=E,I=I,linE=linE,linI=linI,ll=ll)
    end

    event_rates(
        ;β,σ,γ,ω,ψ,pop,S,E,I,R,linE,linI,_...,
    ) = begin
        @assert (I ≥ length(linI) && E ≥ length(linE)) "two: ($I,$(length(linI)),$E,$(length(linI)))"
        alpha = [
            β*S*I/pop,
            β*S*I/pop,
            σ*E,
            σ*E,
            γ*I*indic(I > length(linI)),
            ω*R,
        ]
        pi = [
            (I > 0) ? 1-length(linI)/I : 0,
            (I > 0) ? length(linI)/I : 0,
            (E > 0) ? 1-length(linE)/E : 0,
            (E > 0) ? length(linE)/E : 0,
            1,
            1
        ]
        penalty = ψ*I + γ*I*indic(I ≤ length(linI))
        alpha,pi,penalty
    end

    seir_regular(
        ;t,dt,β,σ,γ,ω,ψ,pop,S,E,I,R,
        ll,linE,linI,_...,
    ) = begin
        tf = t+dt
        if t < tf
            alpha,pi,penalty = event_rates(
                ;β=β,σ=σ,γ=γ,ω=ω,ψ=ψ,pop=pop,S=S,E=E,I=I,R=R,linE=linE,linI=linI,
            )
            @assert (I ≥ length(linI) && E ≥ length(linE)) "three: ($I,$(length(linI)),$E,$(length(linI)))"
            k,s = rcateg(alpha.*pi)
            step::PhyloPOMP.Time = -log(rand())/s
            while (t+step < tf)
                ll -= penalty*step+log(pi[k])
                if k==1
                    S -= 1; E += 1;
                    ll += log(1-length(linE)/E)
                elseif k==2
                    ll += log(length(linI))
                    b = rand(linI); delete!(linI,b); push!(linE,b)
                    S -= 1; E += 1;
                    ll += log(1-length(linI)/I)-log(E)
                elseif k==3
                    E -= 1; I += 1;
                    ll += log(1-length(linI)/I)
                elseif k==4
                    ll += log(length(linE))
                    b = rand(linE); delete!(linE,b); push!(linI,b)
                    E -= 1; I += 1;
                    ll -= log(I)
                elseif k==5
                    I -= 1; R += 1;
                elseif k==6
                    R -= 1; S += 1;
                else
                    @assert false "impossible error!"
                end
                @assert (I ≥ length(linI) && E ≥ length(linE)) "four: ($I,$(length(linI)),$E,$(length(linI)))"
                t += step
                alpha,pi,penalty = event_rates(
                    ;β=β,σ=σ,γ=γ,ω=ω,ψ=ψ,pop=pop,S=S,E=E,I=I,R=R,linE=linE,linI=linI,
                )
                k,s = rcateg(alpha.*pi)
                step = -log(rand())/s
            end
            step = tf - t
            ll -= penalty*step;
        end
        @assert (I ≥ length(linI) && E ≥ length(linE)) "five: ($I,$(length(linI)),$E,$(length(linI)))"
        (S=S,E=E,I=I,R=R,linE=linE,linI=linI,ll=ll)
    end

    seir = function(
        gen::Genealogy;
        β = 4.0, σ = 1.0, γ = 1.0, ω = 1.0, ψ = 0.02,
        pop = 100,
        S0 = 0.9, E0 = 0.0, I0 = 0.02, R0 = 0.08,
        )
        pomp(
            fill((;),length(gen)),
            params = (
                β=Float64(β),σ=Float64(σ),γ=Float64(γ),
                ω=Float64(ω), ψ=Float64(ψ),
                pop=Float64(pop),
                S0=Float64(S0),E0=Float64(E0),
                I0=Float64(I0),R0=Float64(R0),
            ),
            t0 = timezero(gen),
            times = times(gen),
            rinit = function (;S0,E0,I0,R0,pop,_...,)
                m = pop/(S0+E0+I0+R0)
                (
                    S=round(Int64,m*Float64(S0)),
                    E=round(Int64,m*Float64(E0)),
                    I=round(Int64,m*Float64(I0)),
                    R=round(Int64,m*Float64(R0)),
                    node=one(PhyloPOMP.Name),
                    ll=zero(Float64),
                    linE=Set{PhyloPOMP.Name}(),
                    linI=Set{PhyloPOMP.Name}(),
                )
            end,
            rprocess = onestep(
                function(
                    ;S,E,I,R,
                    node,ll,linE,linI,geneal,
                    args...,
                    )
                    linE = copy(linE)
                    linI = copy(linI)
                    ll = zero(Float64)
                    @assert (I ≥ length(linI) && E ≥ length(linE)) "six: $node,($I,$(length(linI)),$E,$(length(linE)))"
                    if E < length(linE)
                        ll += Float64(-Inf)
                        E = length(linE)
                    end
                    if I < length(linI)
                        ll += Float64(-Inf)
                        I = length(linI)
                    end
                    S,E,I,linE,linI,ll = seir_singular(
                        geneal,node;
                        S=S,E=E,I=I,linE=linE,linI=linI,ll=ll,args...,
                    )
                    S,E,I,R,linE,linI,ll = seir_regular(
                        ;S=S,E=E,I=I,R=R,ll=ll,linE=linE,linI=linI,args...,
                    )
                    (S=S,E=E,I=I,R=R,node=node+1,ll=ll,linE=linE,linI=linI)
                end
            ),
            logdmeasure = function (;ll,_...)
                ll
            end,
            userdata=(geneal=gen,),
        )
    end

    p = seir(g)
    @test p isa POMP.PompObject
    @time pf = pfilter(p,Np=100)
    @test isfinite(pf.logLik)
    @test pf isa POMP.PfilterdPompObject

end
