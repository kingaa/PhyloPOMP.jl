module PhyloPOMP

import PartiallyObservedMarkovProcesses as POMP

const Name = UInt64
const Size = UInt64
const Time = Float64

export @demes
include("demes.jl")

export GenealNode, Genealogy
include("genealogy.jl")

include("print.jl")

export parse_newick
include("parse.jl")

export newick
include("newick.jl")

include("fsmarkov.jl")
using .FSMarkov: fsmarkov, transition, generator, statdist
export fsmarkov, transition, generator, statdist

end # module PhyloPOMP
