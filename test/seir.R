library(pomp)
library(phylopomp)
stopifnot(packageVersion("phylopomp") >= "0.19.1")

simulate(
  "SEIR",
  Beta=4,gamma=1,sigma=1,omega=1,psi=0.02,time=50,
  pop=100,S0=0.9,E0=0.0,I0=0.02,R0=0.08,
  ) -> x
x |> plot(points=TRUE)

x |> newick(obscure=TRUE) |> writeLines("seir1.nwk")

readLines("seir1.nwk") |>
  parse_newick(time=50) -> x

tic <- Sys.time()
x |>
  seirs_pomp(
    Beta=4,gamma=1,sigma=1,omega=1,psi=0.02,
    pop=100,S0=0.9,E0=0.0,I0=0.02,R0=0.08
  ) |>
  pfilter(Np=1000) |>
  logLik()
toc <- Sys.time()
toc-tic
