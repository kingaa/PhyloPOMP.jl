module PhyloPOMP

import PartiallyObservedMarkovProcesses as POMP

const Name = UInt64
const Size = UInt64
const Time = Float64
const Prob = Float64

export @demes
include("demes.jl")

export GenealNode, Genealogy
include("genealogy.jl")

export parse_newick
include("parse.jl")

export newick
include("newick.jl")

export fsmarkov
include("fsmarkov.jl")

export guide, relhaz
include("guide.jl")

include("print.jl")

end # module PhyloPOMP
