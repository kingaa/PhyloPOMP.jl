# MERS model

## Naïve filter

```@autodocs
Modules = [PhyloPOMP.NaiveMERS]
Order   = [:module, :function, :macro, :type, :constant]
```

An example constructing the pomp object and running a particle filter.

```@example
using Random: seed!
using PhyloPOMP
using PhyloPOMP.NaiveMERS

seed!(351956486)

p = NaiveMERS.filter_pomp(
    βcc = 2.5, βch = 0.5, βhc = 1.2, βhh = 1.4,
    γc = 1.0, γh = 1.0, χc = 1.0, χh = 1.0,
    Bc = 10.0, Bh = 10.0,
    Sc0 = 0.999, Sh0 = 1.0,
    Ic0 = 0.001, Ih0 = 0.0,
    Nc = 10000.0, Nh = 10000.0,
)
pf = pfilter(p, Np = 1000)
round(logLik(pf),digits=1)
```

## Guided filter

```@autodocs
Modules = [PhyloPOMP.GuidedMERS]
Order   = [:module, :function, :macro, :type, :constant]
```

```@example
using Random: seed!
using PhyloPOMP
using PhyloPOMP.GuidedMERS
using PhyloPOMP.GuidedMERS.Demes: Camel, Human

seed!(351956486)

m = fsmarkov(Camel=>0.5, Human=>0.5, (Camel,Human)=>0.01)
g = parse_newick(GuidedMERS.mers_newick, demes=GuidedMERS.Demes)
p = GuidedMERS.filter_pomp(
	g, m,
	βcc = 2.5, βch = 0.5, βhc = 1.2, βhh = 1.4,
	γc = 1.0, γh = 1.0, χc = 1.0, χh = 1.0,
	Bc = 10.0, Bh = 10.0,
	Sc0 = 0.999, Sh0 = 1.0,
	Ic0 = 0.001, Ih0 = 0.0,
	Nc = 10000.0, Nh = 10000.0,
)
pf = pfilter(p, Np = 1000, trigger=0.2, target=0.8)
round(logLik(pf),digits=1)
```
