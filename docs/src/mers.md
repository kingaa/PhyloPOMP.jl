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

p = NaiveMERS.filter_pomp()
pf = pfilter(p, Np = 1000)
round(logLik(pf),digits=1)
```
