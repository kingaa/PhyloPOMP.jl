using PhyloPOMP
using Test
using Crayons

h1 = crayon"bold blue"
h2 = s -> crayon"!bold light_yellow"("- "*s)

@testset verbose=true "PhyloPOMP.jl" begin

    include("parse.jl")
    include("newick.jl")
    include("cblv.jl")
    include("fsmarkov.jl")
    include("guide.jl")
    include("rcateg.jl")
    include("seir_macro_equivalence.jl")
    include("population_ir_test.jl")
    include("kli_phi_test.jl")
    include("kli_full_transitions_test.jl")
    include("kli_reduce_test.jl")
    include("kli_filter_ir_test.jl")
    include("mgpaudit_test.jl")
    include("kli_decay_test.jl")
    include("kli_proposal_test.jl")
    include("seir_simulate.jl")
    include("kli_seir_compiled_test.jl")
    include("seir_naive.jl")
    include("seir_soft.jl")
    include("seir_guided.jl")
    include("seir_hard.jl")
    include("mers_naive.jl")
    include("mers_simulate.jl")
    include("kli_mers_compiled_test.jl")
    include("mers_soft.jl")
    include("mers_guided.jl")
    include("mers_hard.jl")
    include("kli_properties_test.jl")
    include("kli_kingman_moran_test.jl")
    include("kli_si2r_test.jl")

end
