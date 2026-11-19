public SoftSEIR

"""
    SoftSEIR

An implementation of the phylodynamic filter for an SEIR model, using
a filter guide and "soft" proposals.  That is, when population-process
events occur, the guide can steer those events preferentially onto (or
away from) particular branches, but the overall event rate cannot
increase above that of the underlying population process.
"""
const SoftSEIR = @mgp_soft_filter(
    name=SoftSEIR,
    state=(S,E,I,R),
    demes=(Expos,Infec),
    deme_occupancy=(E,I),
    enns=(nE,nI),
    ells=(ℓE,ℓI),
    rinit = function(; S0, E0, I0, R0, pop, _...)
        barycentric(Int64,(S=S0,E=E0,I=I0,R=R0),pop)
    end,
    init_state = :(S=Int64(0),E=Int64(0),I=Int64(0),R=Int64(0)),
    knowledge = function(v; deme, type, time)
        if type == Sample || type == Node
            demekron!(v, Infec)
            true
        else
            false
        end
    end,
    jumps=(
        trans=(
            rate=β*S*I/pop,
            popmove=begin
                S -= 1
                E += 1
            end,
            moves=(
                reg=(
                    type=swap,
                    from=Infec,
                    to=Expos,
                    idboost=@binomratio(nE,ℓE,1,0),
                    swapboost=@binomratio((nE,nI),(ℓE,ℓI),(1,1),(1,0))
                ),
                fork=(
                    type=fork,
                    from=Infec,
                    to=(Expos,Infec),
                    boost=@binomratio((nE,nI),(ℓE,ℓI),(1,1),(1,1))
                ),
            ),
        ),
        prog=(
            rate=σ*E,
            popmove=begin
                E -= 1
                I += 1
            end,
            moves=(
                swap=(
                    type=swap,
                    from=Expos,
                    to=Infec,
                    idboost=@binomratio(nI,ℓI,1,0),
                    swapboost=@binomratio(nI,ℓI,1,1),
                ),
            ),
        ),
        recovery=(
            rate=γ*I,
            popmove=begin
                I -= 1
                R += 1
            end,
            moves=(
                regular=(
                    type=id,
                    pi=begin
                        @indicator(I>ℓI, 1-ℓI/I)
                    end,
                ),
            ),
        ),
        waning=(
            rate=ω*R,
            popmove=begin
                R -= 1
                S += 1
            end,
            moves=(
                id=(
                    type=id,
                ),
            ),
        ),
        sample_nondestruct=(
            rate=ψ*I,
            moves=(
                inline=(
                    type=chop,
                    signature=(1,1),
                    from=Infec,
                    to=Infec,
                    boost=@binomratio(nI,ℓI,1,1),
                ),
                terminal=(
                    type=chop,
                    signature=(1,0),
                    from=Infec,
                    boost=@binomratio(nI,ℓI,1,0),
                ),
            ),
        ),
        sample_destruct=(
            rate=χ*I,
            popmove=begin
                I -= 1
            end,
            moves=(
                terminal=(
                    type=chop,
                    signature=(1,0),
                    pi=1,
                    from=Infec,
                    boost=1,
                ),
            ),
        ),
    ),
)

SoftSEIR.include("seir_trees.jl")

SoftSEIR.eval(
    :(geneal = [
        parse_newick(seir_trees[1],time=50,demes=DemeSet),
        parse_newick(seir_trees[2],time=400,demes=DemeSet),
    ])
)
