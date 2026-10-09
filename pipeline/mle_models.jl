# Maximum likelihood by mif on the generated filter, in the deep-learning targets' terms, for BD (LBDP), BDEI, BDSS.
# p_sample is known and fixed, as the networks are given it. Each model: its table, the map from the targets to the
# table's rates and back, the founder, and a perturbation for mif that keeps p (μ/χ = (1 − p)/p) and, for BDSS,
# Voznica's constraint λ_ss·λ_nn = λ_sn·λ_ns. Both are linear in the logarithms of the rates, so multiplicative
# perturbations that respect them keep them, and so does mif's average (the geometric mean).

using PhyloPOMP
import PartiallyObservedMarkovProcesses as POMP
using DataFrames, Optim
using Random: AbstractRNG

const MLE_SD = 0.02    # perturbation sd on the log scale, times mif's cooling scale

struct MLEModel
    table::MGPModel
    targets::Vector{String}
    box::Dict{String,Tuple{Float64,Float64}}   # prior ranges of the targets (BDEI's incubation period: 0.2 to 50)
    rates::Function          # (targets vector, p) -> NamedTuple of table parameters
    targets_of::Function     # NamedTuple of rates -> targets vector
    founder::Function        # targets vector -> (x0, graft) for the filter's initial state
    perturb::Function        # mif perturbation
end

pert(scale) = exp(scale * MLE_SD * randn())

const MLE_MODELS = Dict(
    "LBDP" => MLEModel(PhyloPOMP.LBDP, ["R_0", "Infectious_Period"],
        Dict("R_0" => (1.0, 5.0), "Infectious_Period" => (1.0, 10.0)),
        (v, p) -> (γ = 1 / v[2]; (λ = v[1] * γ, μ = γ * (1 - p), ψ = 0.0, χ = γ * p)),
        θ -> (γ = θ.μ + θ.χ; [θ.λ / γ, 1 / γ]),
        v -> ((n = 1,), [1]),
        (scale, lag; λ, μ, χ, _...) -> (f = pert(scale); (λ = λ * pert(scale), μ = μ * f, χ = χ * f))),
    "BDEI" => MLEModel(PhyloPOMP.BDEI, ["R_0", "Infectious_Period", "Incubation_Period"],
        Dict("R_0" => (1.0, 5.0), "Infectious_Period" => (1.0, 10.0), "Incubation_Period" => (0.2, 50.0)),
        (v, p) -> (γ = 1 / v[2]; (σ = 1 / v[3], λ = v[1] * γ, μ = γ * (1 - p), χ = γ * p)),
        θ -> (γ = θ.μ + θ.χ; [θ.λ / γ, 1 / γ, 1 / θ.σ]),
        v -> ((E = 0, I = 1), [0, 1]),
        (scale, lag; σ, λ, μ, χ, _...) -> (f = pert(scale); (σ = σ * pert(scale), λ = λ * pert(scale), μ = μ * f, χ = χ * f))),
    "BDSS" => MLEModel(PhyloPOMP.BDSS, ["R_0", "Infectious_Period", "X_Transmission", "SS_Fraction"],
        Dict("R_0" => (1.0, 5.0), "Infectious_Period" => (1.0, 10.0), "X_Transmission" => (3.0, 10.0),
             "SS_Fraction" => (0.05, 0.2)),
        (v, p) -> begin
            R0, ip, X, f = v
            γ = 1 / ip; B = R0 * γ / (f * X + 1 - f)
            (λ_nn = (1 - f) * B, λ_ns = f * B, λ_sn = X * (1 - f) * B, λ_ss = X * f * B, μ = γ * (1 - p), χ = γ * p)
        end,
        θ -> begin
            γ = θ.μ + θ.χ; B = θ.λ_nn + θ.λ_ns; f = θ.λ_ns / B; X = θ.λ_sn / θ.λ_nn
            [(f * X + 1 - f) * B / γ, 1 / γ, X, f]
        end,
        ## the founder's type is not known from the tree: start from a normal host (the filter's root weights
        ## then decide where the root lineage sits)
        v -> ((N = 1, S = 0), [1, 0]),
        (scale, lag; λ_nn, λ_ns, λ_sn, λ_ss, μ, χ, _...) -> begin
            a, b, c, g = pert(scale), pert(scale), pert(scale), pert(scale)
            (λ_nn = λ_nn * a, λ_ns = λ_ns * b, λ_sn = λ_sn * a * c, λ_ss = λ_ss * b * c, μ = μ * g, χ = χ * g)
        end),
)

"""
    mif_mle(g, mm, p; nstart, Nmif, Np, Np_eval, nreps, maxpop, rng) -> (estimate, loglik, se)

`nstart` mif runs from starts drawn uniformly in the prior box (`profile` runs them in parallel), each end point
evaluated by `nreps` pfilters of `Np_eval` particles; the best end point, in the targets' terms.
"""
function mif_mle(g::Genealogy, mm::MLEModel, p::Real; nstart = 4, Nmif = 40, Np = 500, Np_eval = 1000, nreps = 5,
                 maxpop = 100 * nsample(g), rng::AbstractRNG)
    starts = [[mm.box[t][1] + (mm.box[t][2] - mm.box[t][1]) * rand(rng) for t in mm.targets] for _ in 1:nstart]
    x0, _ = mm.founder(starts[1])
    G = mgp_filter_pomp(g, mm.table; θ = mm.rates(starts[1], p), x0 = x0, maxpop = maxpop)
    design = DataFrame([mm.rates(s, p) for s in starts])
    fit = POMP.profile(G, design; profiled = (), Nmif, Np, Np_eval, nreps, perturbations = mm.perturb,
                       cooling = geometric_cooling(0.5))
    fin = filter(:loglik => isfinite, fit)
    nrow(fin) == 0 && return fill(NaN, length(mm.targets)), -Inf, NaN
    b = fin[argmax(fin.loglik), :]
    θ = NamedTuple{Tuple(keys(mm.rates(starts[1], p)))}(Tuple(b[k] for k in keys(mm.rates(starts[1], p))))
    mm.targets_of(θ), b.loglik, b.se
end

## The founder's deme distribution for the exact likelihood (`mtbd_exact`): generate_trees.jl starts BDSS from a
## superspreader with probability f and the other models from one infectious host.
founder_probs(mm::MLEModel, v) =
    mm.table === PhyloPOMP.BDSS ? [1 - v[4], v[4]] : mm.table === PhyloPOMP.BDEI ? [0.0, 1.0] : [1.0]

## optimizer coordinates: logit for SS_Fraction (a probability), log for the other targets
to_x(t, v) = t == "SS_Fraction" ? log(v / (1 - v)) : log(v)
to_v(t, x) = t == "SS_Fraction" ? 1 / (1 + exp(-x)) : exp(x)

## per host and time unit; inside the prior boxes the largest is a BDSS superspreader's, about 36 (X 10, f 0.05, R0 5, IP 1)
const MTBD_MAXRATE = 100.0

"""
    exact_mle_mtbd(g, mm, p; rtol = 1e-8) -> (estimate, loglik)

Maximum likelihood on the exact likelihood (`mtbd_exact`, not conditioned on the number of samples), p_sample fixed:
Nelder–Mead over the targets (log scale; SS_Fraction on the logit scale) from the centre of the prior box and from
its quarter and three-quarter points. Returns the best end point, in the targets' terms, and its log-likelihood.
Points where a host's total event rate exceeds `MTBD_MAXRATE` count as infinitely bad: the ODE's cost grows with
rate × span, and Nelder–Mead's trial steps otherwise reach incubation periods of 1e-5 (σ ≈ 10⁵, seconds per call).
"""
function exact_mle_mtbd(g::Genealogy, mm::MLEModel, p::Real; rtol = 1e-8)
    ts = mm.targets
    v_of(x) = [to_v(t, xi) for (t, xi) in zip(ts, x)]
    ll(v) = mtbd_exact(g, mm.table; θ = mm.rates(v, p), founder = founder_probs(mm, v), rtol)
    function nll(x)
        v = v_of(x)
        R = mtbd_rates(mm.table, mm.rates(v, p))
        maximum(sum(R.b; dims = 2) .+ sum(R.m; dims = 2) .+ R.d .+ R.ψ .+ R.χ) > MTBD_MAXRATE && return Inf
        l = ll(v)
        isnan(l) ? Inf : -l                                   # a NaN minimum would otherwise win the argmin
    end
    lo = [to_x(t, mm.box[t][1]) for t in ts]
    hi = [to_x(t, mm.box[t][2]) for t in ts]
    fits = [optimize(nll, lo .+ q .* (hi .- lo), NelderMead(), Optim.Options(g_tol = 1e-10, iterations = 10_000))
            for q in (0.5, 0.25, 0.75)]
    best = fits[argmin(Optim.minimum.(fits))]
    v = v_of(Optim.minimizer(best))
    v, ll(v)
end
