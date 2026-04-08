module PhyloPOMP

import PartiallyObservedMarkovProcesses as POMP

include("demes.jl")

export GenealNode, Genealogy
include("genealogy.jl")

export parse_newick
include("parse.jl")

export newick
include("newick.jl")

end # module PhyloPOMP
