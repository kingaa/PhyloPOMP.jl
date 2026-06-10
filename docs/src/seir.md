# SEIR model

## Naïve filter

```@autodocs
Modules = [PhyloPOMP.NaiveSEIR]
Order   = [:module, :function, :macro, :type, :constant]
```

An example constructing the pomp object and running a particle filter.

```@example
using Random: seed!
using PhyloPOMP
using PhyloPOMP.NaiveSEIR

seed!(351956486)

g = parse_newick(seir_trees[1], time = 50.0)
p = seir(g)
pf = pfilter(p, Np = 1000)
round(logLik(pf),digits=1)
```

## Guided filter

```@autodocs
Modules = [PhyloPOMP.GuidedSEIR]
Order   = [:module, :function, :macro, :type, :constant]
```

An example constructing the pomp object and running a particle filter.

```@example
using Random: seed!
using PhyloPOMP
using PhyloPOMP.GuidedSEIR
using PhyloPOMP.GuidedSEIR.Demes: Expos, Infec

seed!(351956486)

g1 = parse_newick(seir_trees[1], time = 50.0)

g2 = guide(
    g1,
    fsmarkov(Expos=>0.1,Infec=>1,(Expos,Infec)=>1),
    seir_convert!
)

p = seir(g2)
pf = pfilter(p, Np = 1000)
round(logLik(pf),digits=1)
```


## Hard-guided filter

```@autodocs
Modules = [PhyloPOMP.HardSEIR]
Order   = [:module, :function, :macro, :type, :constant]
```

An example constructing the pomp object and running a particle filter.

```@example
using Random: seed!
using PhyloPOMP
using PhyloPOMP.HardSEIR
using PhyloPOMP.HardSEIR.Demes: Expos, Infec

seed!(351956486)

g1 = parse_newick(seir_trees[1], time = 50.0)

g2 = guide(
    g1,
    fsmarkov(Expos=>0.1,Infec=>1,(Expos,Infec)=>1),
    seir_convert!
)

p = seir(g2)
pf = pfilter(p, Np = 1000)
round(logLik(pf),digits=1)
```
