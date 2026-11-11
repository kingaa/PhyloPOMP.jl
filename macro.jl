using PhyloPOMP
using PartiallyObservedMarkovProcesses

spec = quote @mgp_soft_filter(
    rinit = function(; S0, E0, I0, R0, pop, _...)
        barycentric(Int64,(S=S0,E=E0,I=I0,R=R0),pop)
    end,
    state=(S,E,I,R),
    demes=(Expos,Infec),
    deme_occupancy=(E,I),
    enns=(nE,nI),
    ells=(ℓE,ℓI),
    live_condition=begin
        nI ≥ ℓI && nE ≥ ℓE
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
                    production=(1,0),
                    saturation=(1,0),

                ),
                fork=(
                    type=fork,
                    from=Infec,
                    to=(Expos,Infec),
                    production=(1,1),
                    saturation=(1,1),
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
                    pi=1,
                    from=Expos,
                    to=Infec,
                    production=(0,1),
                    saturation=(0,1),
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
                    production=(0,0),
                    saturation=(0,0),
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
                    pi=1,
                    type=id,
                    production=(0,0),
                    saturation=(0,0),
                ),
            ),
        ),
        sample_nondestruct=(
            rate=ψ*I,
            moves=(
                inline=(
                    type=chop,
                    signature=(1,1),
                    pi=1,
                    from=Infec,
                    to=Infec,
                    production=(0,1),
                    saturation=(0,1),
                ),
                terminal=(
                    type=chop,
                    signature=(1,0),
                    pi=1,
                    from=Infec,
                    production=(0,1),
                    saturation=(0,0),
                ),
            ),
        ),
        sample_destruct=(
            rate=χ*I,
            moves=(
                terminal=(
                    type=chop,
                    signature=(1,0),
                    pi=1,
                    from=Infec,
                    production=(0,0),
                    saturation=(0,0),
                ),
            ),
        ),
    ),
    )
end

x = eval(spec)

@eval(@macroexpand $spec)
