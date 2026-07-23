#!/usr/bin/env Rscript
# Simulate a forest of SEIR genealogies with phylopomp::runSEIR,
# write Newick + metadata, and plot proof panels.
#
# Usage:
#   Rscript scripts/simulate_seir_forest.R [n_trees] [out_dir] [seed]
#
# Output (under out_dir):
#   forest_manifest.csv   — one row per accepted tree
#   trees/tree_XXXX.nwk   — extended Newick per tree
#   forest_proof.png      — plot grid of the forest
#   PROOF.txt             — human-readable proof of fresh simulation

suppressPackageStartupMessages({
  library(phylopomp)
})

args    <- commandArgs(trailingOnly = TRUE)
n_trees <- if (length(args) >= 1) as.integer(args[1]) else 8L
out_dir <- if (length(args) >= 2) args[2] else
             file.path("output", "seir_forest")
seed    <- if (length(args) >= 3) as.integer(args[3]) else as.integer(Sys.time())

if (!startsWith(out_dir, "/")) {
  # resolve relative to PhyloPOMP.jl root (parent of scripts/)
  script_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) getwd()
  )
  # when invoked as Rscript path/to/script.R, normalize against that path
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  }
  out_dir <- normalizePath(file.path(script_dir, "..", out_dir), mustWork = FALSE)
}

dir.create(file.path(out_dir, "trees"), recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

pkg_ver <- as.character(packageVersion("phylopomp"))
started <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

cat(sprintf(
  "phylopomp %s | seed=%d | n_trees=%d | started=%s\n",
  pkg_ver, seed, n_trees, started
))
stopifnot(package_version(pkg_ver) >= "0.19.0")

# Flu-ish regime + destructive sampling (Voznica-aligned: psi=0)
# Single seed (I0=1/pop) => single-root trees suitable for CBLV.
# Also emit ONE multi-root forest demo so nroot>1 is visible.
Beta  <- 4
sigma <- 1
gamma <- 1
p_sample <- 0.4
chi   <- gamma * p_sample / (1 - p_sample)
omega <- 0
psi   <- 0
pop   <- 1000
S0 <- 0.999; E0 <- 0; I0 <- 0.001; R0 <- 0  # one seed case

min_tips <- 20L
max_tips_single <- 120L
max_tips_multi  <- 400L
max_time <- 40
timeout_tries <- 80L

sim_one <- function(i, multi_root = FALSE) {
  ivp_I0 <- if (multi_root) 0.01 else I0   # ~10 seeds => multi-root forest
  max_tips <- if (multi_root) max_tips_multi else max_tips_single
  # shorter horizon for multi-root so tip counts stay manageable
  t_lo <- if (multi_root) 3 else 8
  t_hi <- if (multi_root) 12 else max_time
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    if (attempt > timeout_tries)
      stop(sprintf("tree %d: exceeded %d attempts", i, timeout_tries))

    time_horizon <- runif(1, t_lo, t_hi)
    x <- tryCatch(
      runSEIR(
        time = time_horizon, t0 = 0,
        Beta = Beta, sigma = sigma, gamma = gamma,
        psi = psi, chi = chi, omega = omega,
        S0 = S0, E0 = E0, I0 = ivp_I0, R0 = R0, pop = pop
      ),
      error = function(e) NULL
    )
    if (is.null(x)) next

    info <- getInfo(x, newick = TRUE, nsample = TRUE, nroot = TRUE, time = TRUE)
    ns <- info$nsample
    nr <- info$nroot
    if (is.null(ns) || is.null(nr) || ns < min_tips || ns > max_tips) next
    if (!multi_root && nr != 1L) next
    if (multi_root && nr < 2L) next

    return(list(
      index = i,
      multi_root = multi_root,
      n_tips = as.integer(ns),
      n_root = as.integer(nr),
      time = info$time,
      attempts = attempt,
      newick = info$newick,
      object = x
    ))
  }
}

rows <- vector("list", n_trees)
objs <- vector("list", n_trees)

# First tree deliberately multi-root (a true genealogical forest)
cat("Simulating tree 1 (multi-root forest demo)...\n")
rows[[1]] <- sim_one(1L, multi_root = TRUE)
objs[[1]] <- rows[[1]]$object
rows[[1]]$object <- NULL
cat(sprintf("  tips=%d roots=%d attempts=%d\n",
            rows[[1]]$n_tips, rows[[1]]$n_root, rows[[1]]$attempts))

# Remaining trees: single-root (CBLV-ready)
for (i in seq.int(2L, n_trees)) {
  cat(sprintf("Simulating tree %d (single-root)...\n", i))
  rows[[i]] <- sim_one(i, multi_root = FALSE)
  objs[[i]] <- rows[[i]]$object
  rows[[i]]$object <- NULL
  cat(sprintf("  tips=%d roots=%d attempts=%d\n",
              rows[[i]]$n_tips, rows[[i]]$n_root, rows[[i]]$attempts))
}

# Write Newick files + SHA1 digests for proof
sha1 <- function(s) digest::digest(s, algo = "sha1", serialize = FALSE)

have_digest <- requireNamespace("digest", quietly = TRUE)
manifest <- data.frame(
  index = integer(),
  multi_root = logical(),
  n_tips = integer(),
  n_root = integer(),
  time = numeric(),
  attempts = integer(),
  nwk_file = character(),
  nwk_sha1 = character(),
  nwk_bytes = integer(),
  stringsAsFactors = FALSE
)

for (i in seq_len(n_trees)) {
  r <- rows[[i]]
  fname <- sprintf("tree_%04d.nwk", i)
  fpath <- file.path(out_dir, "trees", fname)
  writeLines(r$newick, fpath)
  digest_val <- if (have_digest) sha1(r$newick) else NA_character_
  manifest <- rbind(manifest, data.frame(
    index = r$index,
    multi_root = r$multi_root,
    n_tips = r$n_tips,
    n_root = r$n_root,
    time = r$time,
    attempts = r$attempts,
    nwk_file = file.path("trees", fname),
    nwk_sha1 = digest_val,
    nwk_bytes = nchar(r$newick, type = "bytes"),
    stringsAsFactors = FALSE
  ))
}

manifest$Beta <- Beta
manifest$sigma <- sigma
manifest$gamma <- gamma
manifest$psi <- psi
manifest$chi <- chi
manifest$p_sample <- p_sample
manifest$omega <- omega
manifest$pop <- pop
manifest$S0 <- S0
manifest$E0 <- E0
manifest$I0 <- ifelse(manifest$multi_root, 0.01, I0)
manifest$R0 <- R0
manifest$seed <- seed
manifest$phylopomp_version <- pkg_ver
manifest$simulated_at <- started

write.csv(manifest, file.path(out_dir, "forest_manifest.csv"),
          row.names = FALSE, quote = FALSE)

# Plot forest (skip if plotting backends fail)
plot_path <- file.path(out_dir, "forest_proof.png")
ok_plot <- tryCatch({
  png(plot_path, width = 1600, height = 1000, res = 120)
  # phylopomp::plot on gpsim objects; fall back to cowplot if available
  if (requireNamespace("cowplot", quietly = TRUE)) {
    plist <- lapply(objs, function(o) plot(o, points = TRUE))
    print(cowplot::plot_grid(plotlist = plist, ncol = 4))
  } else {
    op <- par(mfrow = c(2, ceiling(n_trees / 2)))
    for (o in objs) plot(o, points = TRUE)
    par(op)
  }
  dev.off()
  TRUE
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  cat("plot failed:", conditionMessage(e), "\n")
  FALSE
})

# PROOF.txt
finished <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
proof <- c(
  "PROOF that these Newick files were freshly simulated with phylopomp::runSEIR",
  "==============================================================================",
  sprintf("phylopomp version : %s  (local install; chi supported)", pkg_ver),
  sprintf("R RNG seed        : %d", seed),
  sprintf("started           : %s", started),
  sprintf("finished          : %s", finished),
  sprintf("n_trees           : %d", n_trees),
  sprintf("params            : Beta=%.3g sigma=%.3g gamma=%.3g psi=%.3g chi=%.3g p_sample=%.3g omega=%.3g pop=%d",
          Beta, sigma, gamma, psi, chi, p_sample, omega, pop),
  "",
  "Forest structure:",
  sprintf("  tree 1 is a MULTI-ROOT genealogical forest (nroot=%d, tips=%d)",
          manifest$n_root[1], manifest$n_tips[1]),
  sprintf("  trees 2..%d are SINGLE-ROOT trees for CBLV (nroot all == 1)", n_trees),
  "",
  "Per-tree tip/root counts (must differ across runs if seed changes):",
  paste(sprintf("  tree %02d: tips=%3d roots=%d attempts=%d sha1=%s",
                manifest$index, manifest$n_tips, manifest$n_root,
                manifest$attempts, substr(manifest$nwk_sha1, 1, 12)),
        collapse = "\n"),
  "",
  sprintf("manifest : %s", file.path(out_dir, "forest_manifest.csv")),
  sprintf("newick   : %s", file.path(out_dir, "trees")),
  sprintf("plot     : %s (%s)", plot_path, if (ok_plot) "ok" else "FAILED"),
  "",
  "This is NOT the canned NaiveSEIR.seir_trees from PhyloPOMP.jl."
)
writeLines(proof, file.path(out_dir, "PROOF.txt"))
cat(paste(proof, collapse = "\n"), "\n")
cat(sprintf("\nWrote forest to %s\n", out_dir))
