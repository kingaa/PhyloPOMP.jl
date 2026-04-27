module PhyloPOMP

import PartiallyObservedMarkovProcesses as POMP
using Reexport: @reexport

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

include("fsmarkov.jl")
@reexport using .FSMarkov: fsmarkov

export guide, relhaz
include("guide.jl")

end # module PhyloPOMP
