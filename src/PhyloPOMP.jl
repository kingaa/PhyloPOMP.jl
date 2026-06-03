"""
    PhyloPOMP

Phylodynamic inference based on POMP models (partially observed Markov processes).

Documentation for **PhyloPOMP.jl** v$(pkgversion(@__MODULE__)).
"""
module PhyloPOMP

const Name = UInt64
const Size = UInt64
const Time = Float64
const Prob = Float64

using Reexport: @reexport
@reexport using EnumX: @enumx

import PartiallyObservedMarkovProcesses as POMP
@reexport using PartiallyObservedMarkovProcesses

export @demes
include("demes.jl")

export GenealNode, Genealogy, times, timezero
include("genealogy.jl")

export parse_newick
include("parse.jl")

export newick
include("newick.jl")

export fsmarkov, generator, forward_action, statdist
include("fsmarkov.jl")

export Guide, guide, relhaz, demekron!
include("guide.jl")

export @marks, rcateg
include("rcateg.jl")

export @indicator
include("indicator.jl")

export Coloring, ell, swap!, chop!, fork!, plant!
include("coloring.jl")

end # module PhyloPOMP
