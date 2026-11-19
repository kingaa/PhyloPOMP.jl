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

import PartiallyObservedMarkovProcesses as POMP
@reexport using PartiallyObservedMarkovProcesses

export @demes
include("demes.jl")

export GenealNode, Genealogy, times, timezero, roots, tips,
    nsample, samples, nodes
include("genealogy.jl")

export parse_newick
include("parse.jl")

export newick
include("newick.jl")

export cblv, parse_cblv
include("cblv.jl")

export Coloring, ell, swap!, chop!, fork!, plant!
include("coloring.jl")

export fsmarkov, generator, forward_action, statdist
include("fsmarkov.jl")

export Guide, GuideNode, guide, relhaz, relhaz_alloc, relhaz!,
    sum_relhaz, demekron!, choose_branch, branch_options_prop,
    branch_options_move
include("guide.jl")

export rcateg
include("rcateg.jl")

export @indicator
include("indicator.jl")

export @binomratio
include("binomratio.jl")

export @mgp_soft_filter, isalive, filter_pomp
include("mgp.jl")

include("examples/Examples.jl")

end # module PhyloPOMP
