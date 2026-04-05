module PhyloPOMP

import PartiallyObservedMarkovProcesses as POMP

export GenealNode, Genealogy
include("genealogy.jl")

export parse_newick
include("parse.jl")

end # module PhyloPOMP
