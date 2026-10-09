## Exact genealogy likelihoods for small models by dynamic programming over individual hosts,
## independent of the KLI algebra and of every filter in PhyloPOMP.jl.
##
## State: the status of every host and, for each lineage of the genealogy alive at time t, the host that
## carries it. Between genealogy nodes the mass evolves by every event that leaves no mark on the tree:
## events that would leave a mark (sampling any host, the end of a lineage without a sample) remove mass. When a host
## that carries a lineage infects another, two hypotheses follow: the lineage stays with the infector or
## passes to the new host. At most one of them is consistent with any complete history, so summing both
## counts each history once. At a node the mass is multiplied by the rate of the event the node records.
##
##   include("scripts/exact_dp.jl")
##   lbdp_dp(g; λ, μ, ψ, χ)              # count space; compare with lbdp_exact
##   seir_dp(g; β, σ, γ, ω, ψ, N, x0)    # labeled hosts, N ≤ 6
##   mers_dp(g; θ, Nc, Nh, Ic0)          # labeled camels and humans, no demography (B_c = B_h = 0)
using PhyloPOMP
using PhyloPOMP: Root, Node, Sample
using SparseArrays, LinearAlgebra

## exp(t A) v for a sparse A (columns are source states): substeps h with h‖A‖₁ ≤ 1, Taylor series in each.
function expmv(A::SparseMatrixCSC, v::Vector{Float64}, t::Float64)
    t <= 0 && return copy(v)
    nrm = opnorm(A, 1)
    nrm == 0 && return copy(v)
    m = max(1, ceil(Int, t * nrm))
    h = t / m
    w = copy(v)
    for _ in 1:m
        term = copy(w); acc = copy(w)
        for k in 1:60
            term = (h / k) * (A * term)
            acc .+= term
            norm(term, 1) <= 1e-17 * norm(acc, 1) && break
        end
        w = acc
    end
    w
end

## ---------------------------------------------------------------------------------------------
## LBDP in count space: state (n, ℓ). Hosts are exchangeable, so the labeled masses collapse exactly:
## births from untracked hosts λ(n-ℓ), from tracked hosts λℓ (stay) + λℓ (move), all to (n+1, ℓ).
function lbdp_dp(g::Genealogy; λ, μ, ψ, χ = 0.0, n0::Integer = 1, nmax::Integer = 80)
    idx(n, l) = n * (nmax + 1) + l + 1          # 0 ≤ l ≤ n ≤ nmax
    D = (nmax + 1)^2
    function generator()
        I = Int[]; J = Int[]; V = Float64[]
        for n in 0:nmax, l in 0:n
            i = idx(n, l)
            push!(I, i); push!(J, i); push!(V, -(λ + μ + ψ + χ) * n)
            n < nmax && (push!(I, idx(n + 1, l)); push!(J, i); push!(V, λ * (n + l)))
            n > l && (push!(I, idx(n - 1, l)); push!(J, i); push!(V, μ * (n - l)))
        end
        sparse(I, J, V, D, D)
    end
    A = generator()
    p = zeros(D); p[idx(n0, 0)] = 1.0
    t = timezero(g)
    for k in eachindex(g)
        nd = g[k]
        p = expmv(A, p, nd.slate - t); t = nd.slate
        q = zeros(D)
        for n in 0:nmax, l in 0:n
            m = p[idx(n, l)]; m == 0 && continue
            if nd.type == Root
                n > l && (q[idx(n, l + 1)] += m * (n - l))
            elseif nd.type == Node
                n < nmax && l >= 1 && (q[idx(n + 1, l + 1)] += m * 2λ)
            elseif isempty(nd.children)                       # tip sample
                l >= 1 && (q[idx(n, l - 1)] += m * ψ)
                l >= 1 && n >= 1 && (q[idx(n - 1, l - 1)] += m * χ)
            else                                              # inline sample
                l >= 1 && (q[idx(n, l)] += m * ψ)
            end
        end
        p = q
    end
    p = expmv(A, p, g.time - t)
    log(sum(p))
end

## ---------------------------------------------------------------------------------------------
## SEIR with N labeled hosts. Status 0 = S, 1 = E, 2 = I, 3 = R. Infection: each (I host, S host) pair at
## rate β/N; progression σ per E host; recovery γ per I host; waning ω per R host; sampling ψ per I host
## (the host stays infectious). `drop_stay = true` switches off the "lineage stays with the infector"
## hypothesis whenever every infectious host carries a lineage: the expectation of a filter that never
## proposes that outcome.
function seir_dp(g::Genealogy; β, σ, γ, ω, ψ, N::Integer, x0::NamedTuple, drop_stay::Bool = false)
    st0 = Int8[fill(0, x0.S); fill(1, x0.E); fill(2, x0.I); fill(3, x0.R)]
    length(st0) == N || error("x0 must sum to N")
    lins = Int[]                                   # lineage ids alive, in a fixed order
    dist = Dict{Any,Float64}((Tuple(st0), ()) => 1.0)

    ## transitions out of a state that leave no mark on the tree
    function moves(st, hs)
        out = Tuple{Any,Float64}[]
        carrier = Dict(Int(h) => j for (j, h) in enumerate(hs))
        nI = count(==(2), st)
        alltracked = nI > 0 && count(h -> st[h] == 2, hs) == nI
        exitrate = 0.0
        for h in 1:N
            s = st[h]
            if s == 2
                for k in 1:N
                    st[k] == 0 || continue
                    exitrate += β / N
                    st2 = Base.setindex(st, Int8(1), k)
                    if haskey(carrier, h)
                        (drop_stay && alltracked) || push!(out, ((st2, hs), β / N))       # stays
                        push!(out, ((st2, Base.setindex(hs, Int8(k), carrier[h])), β / N)) # moves
                    else
                        push!(out, ((st2, hs), β / N))
                    end
                end
                exitrate += γ + ψ
                haskey(carrier, h) || push!(out, ((Base.setindex(st, Int8(3), h), hs), γ))
            elseif s == 1
                exitrate += σ
                push!(out, ((Base.setindex(st, Int8(2), h), hs), σ))
            elseif s == 3
                exitrate += ω
                push!(out, ((Base.setindex(st, Int8(0), h), hs), ω))
            end
        end
        out, exitrate
    end

    function propagate!(dist, dt)
        dt <= 0 && return dist
        keys0 = collect(keys(dist))
        index = Dict{Any,Int}(k => i for (i, k) in enumerate(keys0))
        states = copy(keys0)
        I = Int[]; J = Int[]; V = Float64[]
        i = 1
        while i <= length(states)
            (st, hs) = states[i]
            out, ex = moves(st, hs)
            push!(I, i); push!(J, i); push!(V, -ex)
            for (dest, r) in out
                j = get!(index, dest) do
                    push!(states, dest); length(states)
                end
                push!(I, j); push!(J, i); push!(V, r)
            end
            i += 1
        end
        A = sparse(I, J, V, length(states), length(states))
        p = zeros(length(states)); for (k, m) in dist; p[index[k]] = m; end
        p = expmv(A, p, dt)
        empty!(dist)
        for (k, m) in zip(states, p); m > 0 && (dist[k] = m); end
        dist
    end

    t = timezero(g)
    for k in eachindex(g)
        nd = g[k]
        propagate!(dist, nd.slate - t); t = nd.slate
        new = Dict{Any,Float64}()
        add!(key, m) = (new[key] = get(new, key, 0.0) + m)
        if nd.type == Root
            push!(lins, nd.lineage)
            for ((st, hs), m) in dist, h in 1:N
                (st[h] in (1, 2) && !(Int8(h) in hs)) || continue
                add!((st, (hs..., Int8(h))), m)
            end
        else
            j = findfirst(==(nd.lineage), lins)
            j === nothing && error("lineage $(nd.lineage) not alive at node $k")
            if nd.type == Node
                c1, c2 = (g[c].lineage for c in nd.children)
                newlins = [deleteat!(copy(lins), j); c1; c2]
                for ((st, hs), m) in dist
                    h = hs[j]; st[h] == 2 || continue
                    rest = Tuple(hs[i] for i in eachindex(hs) if i != j)
                    for kk in 1:N
                        st[kk] == 0 || continue
                        st2 = Base.setindex(st, Int8(1), kk)
                        add!((st2, (rest..., h, Int8(kk))), m * β / N)   # c1 in infector, c2 in new host
                        add!((st2, (rest..., Int8(kk), h)), m * β / N)   # c1 in new host, c2 in infector
                    end
                end
                lins = newlins
            elseif isempty(nd.children)                                    # tip sample
                for ((st, hs), m) in dist
                    st[hs[j]] == 2 || continue
                    add!((st, Tuple(hs[i] for i in eachindex(hs) if i != j)), m * ψ)
                end
                deleteat!(lins, j)
            else                                                           # inline sample
                c = g[only(nd.children)].lineage
                for ((st, hs), m) in dist
                    st[hs[j]] == 2 || continue
                    rest = Tuple(hs[i] for i in eachindex(hs) if i != j)
                    add!((st, (rest..., hs[j])), m * ψ)
                end
                lins = [deleteat!(copy(lins), j); c]
            end
        end
        dist = new
    end
    propagate!(dist, g.time - t)
    log(sum(values(dist)))
end

## ---------------------------------------------------------------------------------------------
## Mass transport between nodes for any labeled-state model: `moves(state) -> (transitions, exit rate)`.
function _propagate(dist::Dict, dt::Float64, moves)
    dt <= 0 && return dist
    states = collect(keys(dist))
    index = Dict{Any,Int}(k => i for (i, k) in enumerate(states))
    I = Int[]; J = Int[]; V = Float64[]
    i = 1
    while i <= length(states)
        out, ex = moves(states[i])
        push!(I, i); push!(J, i); push!(V, -ex)
        for (dest, r) in out
            j = get!(index, dest) do
                push!(states, dest); length(states)
            end
            push!(I, j); push!(J, i); push!(V, r)
        end
        i += 1
    end
    A = sparse(I, J, V, length(states), length(states))
    p = zeros(length(states)); for (k, m) in dist; p[index[k]] = m; end
    p = expmv(A, p, dt)
    Dict{Any,Float64}(k => m for (k, m) in zip(states, p) if m > 0)
end

## MERS with Nc labeled camels (hosts 1:Nc) and Nh labeled humans, no births or deaths of susceptibles.
## Status 0 = S, 2 = I, 4 = gone (removed or sampled). Infection by an I host h of an S host k at rate
## β_cc/N_c (camel→camel), β_hc/N_c (camel→human), β_hh/N_h (human→human), β_ch/N_h (human→camel);
## removal γ and sampling χ (which removes the host) per I host of each species. Sample demes are observed:
## a tip of deme Camel (Human) must be a camel (human) host. `drop_stay = true` removes the "stays"
## hypothesis of a cross-species infection when every infected host of the infector's species is tracked.
function mers_dp(g::Genealogy; θ, Nc::Integer, Nh::Integer, Ic0::Integer = 1, drop_stay::Bool = false)
    N = Nc + Nh
    camel(h) = h <= Nc
    rate(h, k) = camel(h) ? (camel(k) ? θ.β_cc : θ.β_hc) / θ.N_c : (camel(k) ? θ.β_ch : θ.β_hh) / θ.N_h
    st0 = Int8[i <= Ic0 ? 2 : 0 for i in 1:N]
    lins = Int[]
    dist = Dict{Any,Float64}((Tuple(st0), ()) => 1.0)
    function moves(state)
        st, hs = state
        out = Tuple{Any,Float64}[]
        carrier = Dict(Int(h) => j for (j, h) in enumerate(hs))
        alltracked(sp) = (nI = count(h -> st[h] == 2 && camel(h) == sp, 1:N);
                          nI > 0 && count(h -> st[h] == 2 && camel(h) == sp, hs) == nI)
        ex = 0.0
        for h in 1:N
            st[h] == 2 || continue
            for k in 1:N
                st[k] == 0 || continue
                r = rate(h, k); ex += r
                st2 = Base.setindex(st, Int8(2), k)
                if haskey(carrier, h)
                    cross = camel(h) != camel(k)
                    (drop_stay && cross && alltracked(camel(h))) || push!(out, ((st2, hs), r))
                    push!(out, ((st2, Base.setindex(hs, Int8(k), carrier[h])), r))
                else
                    push!(out, ((st2, hs), r))
                end
            end
            γ, χ = camel(h) ? (θ.γ_c, θ.χ_c) : (θ.γ_h, θ.χ_h)
            ex += γ + χ
            haskey(carrier, h) || push!(out, ((Base.setindex(st, Int8(4), h), hs), γ))
        end
        out, ex
    end
    t = timezero(g)
    for k in eachindex(g)
        nd = g[k]
        dist = _propagate(dist, nd.slate - t, moves); t = nd.slate
        new = Dict{Any,Float64}()
        add!(key, m) = (new[key] = get(new, key, 0.0) + m)
        if nd.type == Root
            push!(lins, nd.lineage)
            for ((st, hs), m) in dist, h in 1:N
                (st[h] == 2 && !(Int8(h) in hs)) && add!((st, (hs..., Int8(h))), m)
            end
        else
            j = findfirst(==(nd.lineage), lins)
            if nd.type == Node
                c1, c2 = (g[c].lineage for c in nd.children)
                for ((st, hs), m) in dist
                    h = hs[j]; st[h] == 2 || continue
                    rest = Tuple(hs[i] for i in eachindex(hs) if i != j)
                    for kk in 1:N
                        st[kk] == 0 || continue
                        st2 = Base.setindex(st, Int8(2), kk); r = rate(h, kk)
                        add!((st2, (rest..., h, Int8(kk))), m * r)
                        add!((st2, (rest..., Int8(kk), h)), m * r)
                    end
                end
                lins = [deleteat!(copy(lins), j); c1; c2]
            else                                                        # tip sample, removes the host
                isempty(nd.children) || error("mers_dp: sampled ancestors are not supported (sampling removes the host)")
                iscamel = string(nd.deme) == "Camel"
                for ((st, hs), m) in dist
                    h = hs[j]
                    (st[h] == 2 && camel(h) == iscamel) || continue
                    add!((Base.setindex(st, Int8(4), h), Tuple(hs[i] for i in eachindex(hs) if i != j)),
                         m * (iscamel ? θ.χ_c : θ.χ_h))
                end
                deleteat!(lins, j)
            end
        end
        dist = new
    end
    dist = _propagate(dist, g.time - t, moves)
    log(sum(values(dist)))
end
