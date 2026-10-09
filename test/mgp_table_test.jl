"""
`filter_table`: one row per event, rate expressions from `@mgp`, filter weights checked against the IR.
"""
module MgpTableTest

import ..Main: h1, h2

@info h1("filter_table")

using Test
using PhyloPOMP

@testset "filter_table" begin
    for M in (PhyloPOMP.SEIR, PhyloPOMP.MERS, SIR, SI2R, LBDP, BDEI, BDSS, PhyloPOMP.MTBD)
        md = filter_table(M)                       # throws if a formula disagrees with the IR
        @test count(==('\n'), md) == length(M.events) + 2
        tex = filter_table(M; format = :latex)
        @test startswith(tex, "\\begin{tabular}") && occursin("\\end{tabular}", tex)
        @test !occursin('β', tex) && !occursin('ℓ', tex)
    end
    md = filter_table(PhyloPOMP.MERS)
    @test occursin("| transmission_hc | BIRTH | (β_hc * S_h * I_c) / N_c |", md)
    @test occursin("no lineage moves: 1 − ℓ_I_h/n′_I_h", md)
    @test occursin("1 − C(ℓ_I_c,2)/C(n′_I_c,2)", md)
    @test occursin("| progression | MIGRATION |", filter_table(PhyloPOMP.SEIR))
    @test occursin("(closure)", filter_table(PhyloPOMP.SEIR_REFERENCE))   # hand-built events keep no expression
    @test_throws ArgumentError filter_table(SIR; format = :html)
end

end # module MgpTableTest
