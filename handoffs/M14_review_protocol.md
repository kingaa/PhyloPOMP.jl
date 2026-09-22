# Review Protocol — M14: Upstream Guided-Filter Sync

## 0. Purpose and how to use this

An independent agent (or human) should use this to verify the work described
in `handoffs/M14_upstream_guided_sync.md`, without trusting any claim in that
file or in the session's chat transcript at face value. Every task below
names the exact file/line, gives a command to run, states the expected
result, and gives explicit PASS/FAIL/CONCERN criteria. Where a check is
stochastic (Monte Carlo logLik estimates), a tolerance is given — do not
treat a near-miss as automatic failure, but do not wave through anything
outside the stated bound either.

Work through tasks in order. Do not skip Task 8 (git hygiene) even if
everything else passes — scope creep is a failure mode of this kind of
multi-agent session independent of correctness.

**Report format:** for each subtask, record PASS / FAIL / CONCERN, the
command run, the actual output (not paraphrased), and a one-line verdict.
Do not summarize away a FAIL or CONCERN by rolling it into an overall PASS.

## 1. Directory map

```
/home/inferno/Desktop/Projects/PhyloPOMP.jl          our repo, branch atpabuser-devel — THE THING BEING REVIEWED
/home/inferno/Desktop/Projects/PhyloPOMP.jl_aarons   read-only reference checkout of upstream (kingaa/PhyloPOMP.jl), branch mers_guide
```

Inside `PhyloPOMP.jl`:
```
src/                          core package (guide.jl, coloring.jl, genealogy.jl, parse.jl, cblv.jl, ...)
src/examples/                 filter implementations (mers_guided.jl, mers_soft.jl, mers_hard.jl,
                               mers_naive.jl, mers_funs.jl, seir_*.jl, mgp*.jl)
test/                         test suite, run via test/runtests.jl or individually
docs/                         documentation, including docs/compiler/*.tex (assessment notes)
handoffs/                     milestone handoff docs (this file's sibling: M14_upstream_guided_sync.md)
```

Nothing in this session was committed. `git status --short` at the repo root
is the ground truth for what changed — cross-check it against §9 below rather
than trusting any file list in prose.

## 2. Task 1 — Dependency and environment sync

### 2.1 Project.toml compat bump
```
git diff Project.toml
```
Expect: `PartiallyObservedMarkovProcesses` compat changed from `"0.8"` to
`"0.10"`, the `[sources]` devel-branch pin removed, `Distributions` added to
`[deps]`, `BenchmarkTools`/`Crayons`/`Tally` removed from `[deps]` (they
should appear only in `[extras]` or `test/Project.toml`), version bumped to
`0.0.13-3`.

**Verify independently, don't just read the diff:**
```
julia --project -e 'using Pkg; Pkg.status()'
```
PASS iff `PartiallyObservedMarkovProcesses v0.10.0` resolves from the
General registry (not a `repo-url`/`repo-rev` devel pin). Check
`Manifest.toml` directly if uncertain:
```
grep -A6 '\[\[deps.PartiallyObservedMarkovProcesses\]\]' Manifest.toml
```
FAIL if `registries = "General"` is absent or a `repo-rev` key is present.

### 2.2 test/Project.toml
```
git diff test/Project.toml
```
Expect `Distributions` added (the guided-MERS test imports `LogNormal`).
PASS iff `julia --project=test -e 'using Distributions'` succeeds with no
`Pkg.add` needed.

### 2.3 CI workflow
```
git diff .github/workflows/CI.yml
diff .github/workflows/CI.yml ../PhyloPOMP.jl_aarons/.github/workflows/CI.yml
```
PASS iff the second diff is empty (byte-identical to upstream) or the only
difference is something explicitly justified in the handoff.

### 2.4 guide.jl and docstring sync
```
diff src/guide.jl ../PhyloPOMP.jl_aarons/src/guide.jl
diff src/genealogy.jl ../PhyloPOMP.jl_aarons/src/genealogy.jl
diff src/parse.jl ../PhyloPOMP.jl_aarons/src/parse.jl
diff src/cblv.jl ../PhyloPOMP.jl_aarons/src/cblv.jl
diff src/coloring.jl ../PhyloPOMP.jl_aarons/src/coloring.jl
```
PASS iff all five diffs are empty. These files are claimed to be byte-for-byte
synced to upstream — any remaining diff is a FAIL, not a rounding error.

## 3. Task 2 — The two `mers_guided.jl` bug fixes (mathematical correctness)

This is the highest-stakes claim in the handoff: that upstream's rewrite of
`src/examples/mers_guided.jl` (their commit `25039d0`) has two likelihood
errors, both patched here. Do not accept the session's own numerical checks
as sufficient — re-derive independently.

### 3.1 Bug 1 — within-host fork factor (lines 187, 219)
```
sed -n '165,225p' src/examples/mers_guided.jl
```
Current code charges `ll -= log(I_c*(I_c-1)/2)` (line 187, camel-camel fork)
and `ll -= log(I_h*(I_h-1)/2)` (line 219, human-human fork).

**Independent derivation:** at a fork where a Camel individual with `I_c`
total camel infecteds (all indistinguishable to the observer) gives rise to
two tracked child lineages, the number of ways to choose an ordered pair of
distinct individuals from `I_c` is `I_c*(I_c-1)`; the number of *unordered*
pairs (the actual sample space for "which two are the observed lineages") is
`I_c*(I_c-1)/2`. The probability the fork is exactly the observed pair is
`1/binomial(I_c,2) = 2/(I_c*(I_c-1))`, i.e. `-log(I_c*(I_c-1)/2)`. If you
derive a different answer, STOP and flag as CONCERN before proceeding — do
not just defer to the code.

**Cross-check against the compiler IR** (independent code path, not written
by the same session's fix):
```
julia --project -e '
using PhyloPOMP
using PhyloPOMP: full_transitions, ForkTransition
tcc = PhyloPOMP.MERS.events[findfirst(e->e.name==:transmission_cc, PhyloPOMP.MERS.events)]
for Ic in (2,3,5,8)
    ts = full_transitions(tcc, [2,0], [Ic,0])
    fk = filter(t->t isa ForkTransition, ts)
    println("I_c=$Ic  phi=", isempty(fk) ? "absent" : fk[1].phi, "  expected=", 2//(Ic*(Ic-1)))
end'
```
PASS iff `phi` equals `expected` at every `I_c` tested.

**Cross-check against the independent hand-coded oracle** (`mers_naive.jl`,
never touched by this session):
```
grep -n 'log(lambda_cc) - log(Ic \* (Ic - 1) / 2)\|log(lambda_hh) - log(Ih \* (Ih - 1) / 2)' src/examples/mers_naive.jl
```
PASS iff both lines are present with the `/2`.

### 3.2 Bug 2 — missing no-fork factor on regular within-host births (lines 267–288)
```
sed -n '260,290p' src/examples/mers_guided.jl
```
Current code: `regular_transmission_cc!`/`regular_transmission_hh!` return
`log(1-ellC*(ellC-1)/I_c/(I_c-1))` (line 275) / `log(1-ellH*(ellH-1)/I_h/(I_h-1))`
(line 286) rather than a zero log-weight.

**Independent derivation:** at a birth in *regular* time (i.e. one that does
not correspond to an observed coalescent node), the event must be weighted by
the probability that the birth did *not* involve two of the `ellC` currently
tracked lineages — otherwise it would have produced an observable branch
point that isn't in the data. That probability is
`1 - binomial(ellC,2)/binomial(I_c,2) = 1 - ellC*(ellC-1)/(I_c*(I_c-1))`.
Confirm this matches the Chu–Vandermonde identity check already in the test
suite:
```
grep -n -B3 -A15 'weighted aggregate\|TCC/THH' test/kli_mers_compiled_test.jl | head -60
```
PASS iff that testset's own algebra (`Φ_identity + ℓ·Φ_inline = 1 - C(ℓ,2)·Φ_fork`)
is logically the same statement as the code fix. Run it:
```
julia --project=test -e '
using Test, Crayons
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
include("test/kli_mers_compiled_test.jl")'
```
PASS iff all its assertions pass (this testset is independent of
`mers_guided.jl` — it tests the compiler IR against `mers_naive.jl`, not
against the fix under review — so it corroborates the *formula*, not the
*specific edit*).

### 3.3 Independent numerical reproduction of both bugs together
Do not reuse the session's scratch scripts (they no longer exist — this must
be written fresh). Minimal reproduction:
```
julia --project -e '
using PhyloPOMP, PhyloPOMP.GuidedMERS
using Random: seed!
nwk = "((([&&PhyloPOMP deme=camel]1:1.0,[&&PhyloPOMP deme=camel]2:1.0):0.5,[&&PhyloPOMP deme=camel]3:1.2):0.3,[&&PhyloPOMP deme=camel]4:1.9);"
g = parse_newick(nwk, demes=GuidedMERS.Demes)
m = fsmarkov(GuidedMERS.Camel=>0.5, GuidedMERS.Human=>0.5, (GuidedMERS.Camel,GuidedMERS.Human)=>0.01)
p = GuidedMERS.filter_pomp(g, m; β_cc=3.0, β_ch=0.0, β_hc=0.0, β_hh=3.0,
    γ_c=1.0, γ_h=1.0, χ_c=1.0, χ_h=1.0, B_c=0.0, B_h=0.0,
    S_c0=1.0, S_h0=1.0, I_c0=0.5, I_h0=0.0, N_c=60, N_h=60)
seed!(99)
ll = [logLik(pfilter(p, Np=4000)) for _ in 1:20]
est, se = logmeanexp(ll, se=true)
println("fixed kernel: ", round(est,digits=3), " ± ", round(se,digits=3))'
```
This tree has 3 internal nodes, all forced camel-camel by `β_hc=β_ch=0`.
Manually revert both fixes in a scratch copy (change line 187/219 to drop the
`/2`, and lines 275/286 to `zero(Prob)`), rerun, and confirm the *fixed*
estimate is higher (less negative) by roughly `3*log(2) ≈ 2.08` from
reverting bug 1 alone, and substantially higher again from also reverting
bug 2 on a tree with long coexisting branches (bug 2 needs `ellC ≥ 2` to
bite, so this 4-tip example alone under-tests it — see §3.1's derivation
instead for that half, or build a tree with 4+ simultaneous camel lineages
before any coalescence).

CONCERN (not necessarily FAIL) if you cannot reproduce a consistent
directional effect — re-read the derivation in §3.1/§3.2 before concluding
the fix itself is wrong, since a sign error in your own revert is equally
possible.

## 4. Task 3 — Soft/Hard MERS kernel refactor

### 4.1 Structural parity check
```
diff <(grep -oE '^[a-zA-Z_!]+\(' src/examples/mers_funs.jl | sort -u) \
     <(grep -oE '^[a-zA-Z_!]+\(' src/examples/seir_funs.jl | sort -u)
```
This won't match 1:1 (MERS has more compartments) but should show the same
*shape* of names: `knowledge!`, `check`, rate pieces, `singular_*!`,
`*_rinit`, `filter_pomp`. FAIL if `mers_funs.jl` is missing any of
`singular_root!`, `terminal_sample!`, `singular_branch!`, `singular_part!`.

### 4.2 Verbatim-copy claim
The handoff claims the singular part of `mers_funs.jl` is copied verbatim
from the (fixed) `mers_guided.jl`. Verify:
```
diff <(sed -n '/^deme_occupancy/,/^end$/p' src/examples/mers_funs.jl) \
     <(sed -n '/^deme_occupancy/,/^end$/p' src/examples/mers_guided.jl)
```
Run the equivalent diff for `singular_root!`, `terminal_sample!`,
`singular_sample!`, `singular_branch!`, `singular_part!`, `mers_rinit`. Minor
whitespace/comment differences are fine; any *logic* difference is a FAIL —
the entire point of the refactor is that all three kernels share this code
exactly.

### 4.3 `event_rates!` on*/off* convention
```
sed -n '/^transmission!/,/^end$/p' src/examples/mers_funs.jl
```
Confirm Soft's `regular_part!` in `mers_soft.jl` passes `onC=ellC, onH=ellH`
and Hard's passes `onC=sum_relhaz(...), offC=I_c-ellC` (grep both files for
`onC=`). FAIL if either kernel passes something structurally different from
what `mers_soft.jl`'s own header comment (lines ~35-46) says it does — a
comment that doesn't match the code below it is itself a defect to report.

### 4.4 Removal-accounting equivalence proof
`mers_soft.jl` carries an inline comment proving the new `removal!`
convention (`pi=1-ell/I`, generic `-log(pi[k])` charge) is value-identical to
the old hand-charged convention. Re-derive this yourself on paper (it's
three lines of algebra: alpha, charge, decay in both conventions) rather than
trusting the comment. FAIL if the two conventions disagree in either the
`I > ell` or `I ≤ ell` case.

### 4.5 Regression against pre-refactor baseline
```
git show 93ebd17~1:src/examples/mers_soft.jl > /tmp/old_mers_soft.jl
```
(That commit predates the whole upstream sync — the file used the *old*
`mers_funs.jl` and old defaults.) Diff its `common` test parameters against
current `test/mers_soft.jl`'s; note that defaults changed (old had
`χ_h=β_hc=β_ch=0`, new shares `GuidedMERS`'s non-degenerate defaults), so a
literal old-vs-new logLik comparison at *default* parameters is not
meaningful — only a like-for-like comparison (explicit, matched kwargs) is.
Confirm the handoff's regression table used explicit matched kwargs, not
bare defaults, before trusting it.

### 4.6 Three-kernel agreement — reproduce independently
```
julia --project=test -e '
using Test, Crayons
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
include("test/mers_soft.jl")'
```
PASS iff the "soft, hard and guided agree on the likelihood" testset passes
(all pairwise z < 3). Then independently, at higher Np (10000, ≥20
replicates per kernel, fresh seeds not reused from the test file), confirm
the three-kernel — plus the *pre-rewrite* `GuidedMERS` from
`git show 93ebd17~1:src/examples/mers_guided.jl` as a fourth, older-vintage
estimator — all agree pairwise within roughly 1 SE. A looser agreement
(z up to 3) is not automatically a FAIL given the fixture's documented heavy
tails (see §4.7), but should be re-run at higher Np before dismissing.

### 4.7 Heavy-tail claim
The handoff and test comments claim `logmeanexp` is heavy-tailed on the
3-tip `small_tree` fixture, and that a *larger* replicate count can paradoxically
give a *worse* worst-case z-score. Spot-check this claim by running the
agreement block at R=10, R=15, R=20 with several seeds and confirming the
qualitative pattern (non-monotonic reliability in R) rather than assuming
it's a rationalization for a marginal result.

## 5. Task 4 — Test suite changes

### 5.1 `test/mers_guided.jl` heavy-gating
```
sed -n '1,65p' test/mers_guided.jl
```
Confirm: (a) the file parses (`julia --project=test -e 'Meta.parseall(read("test/mers_guided.jl",String)); println("OK")'`);
(b) the `mif` block is inside `if heavy ... end`, not gated by an
out-of-function `return` (a `return` at `@testset` top level inside a module
is itself a lowering error — confirm this was actually the bug, not assumed);
(c) `@test pf isa POMP.PfilterdPompObject` runs unconditionally;
(d) `RUN_HEAVY_TESTS=no` genuinely skips the `mif` block in under 15 seconds:
```
cd test && RUN_HEAVY_TESTS=no time julia --project=. -e '
using Test, Crayons
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
include("mers_guided.jl")'
```

### 5.2 `test/mers_simulate.jl` / `test/seir_simulate.jl` additions
```
git diff test/mers_simulate.jl test/seir_simulate.jl
```
Confirm the new testsets actually call `GuidedMERS.check`/`GuidedSEIR.check`
(not just construct the filter) and assert `isfinite(logLik(...))`, not just
`isa PfilterdPompObject` (a type check alone would pass even at `-Inf`).
Re-run both files fresh, at least 3 times with different `ENV["JULIA_SEED"]`
or by editing the seed literal, to rule out a lucky one-off:
```
cd test && for i in 1 2 3; do RUN_HEAVY_TESTS=no julia --project=. -e '
using Test, Crayons, Random
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
include("mers_simulate.jl"); include("seir_simulate.jl")' 2>&1 | grep -E "Fail|Error|Pass.*Total"; done
```

### 5.3 Full suite
```
cd test && RUN_HEAVY_TESTS=no julia --project -e 'using Pkg; Pkg.test()'
```
PASS iff this exits 0 with "tests passed" and no errors. Count files
actually included:
```
grep -c 'include(' test/runtests.jl
```
The handoff claims "30 files" — verify this number directly rather than
trusting it copied correctly.

## 6. Task 5 — Documentation accuracy

### 6.1 `docs/compiler/guided_kernel_assessment.tex` addendum
Read the "Addendum" section added at the end. For each numbered claim in it,
cross-check against the corresponding code/test evidence in Tasks 1–4 above.
FAIL any claim that isn't independently verifiable by a task above — a
document that asserts something no test or diff supports is worse than no
document.

Confirm the PDF was actually rebuilt from the current `.tex` (not stale):
```
cd docs/compiler && pdflatex -interaction=nonstopmode guided_kernel_assessment.tex >/tmp/relatex.log 2>&1
grep -i 'error' /tmp/relatex.log
diff <(md5sum guided_kernel_assessment.pdf) <(pdftotext guided_kernel_assessment.pdf - | md5sum) # sanity only, not a real check
```
Simplest real check: `pdftotext guided_kernel_assessment.pdf -` and grep for
the addendum's key phrases (e.g. "Addendum") to confirm they're actually in
the rendered PDF, not just the source.

### 6.2 `src/examples/mers_filter_suite.tex` changelog
```
git diff src/examples/mers_filter_suite.tex
```
Confirm the new changelog paragraph's numbers match §3.3's reproduction, and
that it does NOT claim the bugs were "already reported to Aaron" (an earlier
draft in this session had that wording; it should say the fixes are *to be*
relayed, not that they *have been*, unless you have independent evidence
Aaron was actually told).

### 6.3 README.md
```
diff README.md ../PhyloPOMP.jl_aarons/README.md
```
PASS iff identical or near-identical (badge/section-order sync only, no
content divergence introduced).

### 6.4 `mers_guided_part1.md`
Confirm the prepended note accurately says this file describes a reverted,
superseded refactor, and does not contradict `handoffs/M14_upstream_guided_sync.md`.

## 7. Task 6 — Independent reproduction of the Aaron's-repo claims

These are strong, falsifiable claims made mid-session. Re-derive them from
scratch — do not just re-read the session's own commands.

### 7.1 "test/mers_guided.jl is missing a closing `end`"
```
cd ../PhyloPOMP.jl_aarons && git fetch origin mers_guide && git log -1 --format='%H' origin/mers_guide
grep -n '^module\|^end' test/mers_guided.jl
julia -e 'Meta.parseall(read("test/mers_guided.jl",String))' 2>&1 | head -5
```
PASS (claim confirmed) iff the parse throws `ParseError: ... Expected `end``.
Record the exact commit hash you tested against — if upstream has moved
since this session, redo this check and note the new commit.

### 7.2 "Distributions is missing from Aaron's test/Project.toml"
```
grep -n 'Distributions' test/Project.toml   # should be absent
grep -n 'using Distributions' test/mers_guided.jl  # should be present
```
PASS (claim confirmed) iff the first grep finds nothing and the second finds
line 10 (`using Distributions: LogNormal`).

### 7.3 "test/mers_guided.jl was never included in test/runtests.jl, at any relevant commit"
```
git show 44edaf6:test/runtests.jl | grep -n mers
git show 25039d0:test/runtests.jl | grep -n mers
grep -n mers test/runtests.jl   # current HEAD of mers_guide
```
PASS (claim confirmed) iff none of the three show an `include("mers_guided.jl")` line.

### 7.4 "CI never actually ran it — confirmed via coverage log, not just absence from runtests.jl"
```
gh run list --repo kingaa/PhyloPOMP.jl --branch mers_guide --limit 5
gh api repos/kingaa/PhyloPOMP.jl/actions/jobs/<job-id>/logs | grep -i "mers_guided.jl\|GuidedMERS\|MERS model with guided"
```
(Use a job ID from the `gh run list` output — pick the most recent
successful run on `mers_guide`.) PASS iff the log shows
`Coverage file(s) for src/examples/mers_guided.jl do not exist` and no
"MERS model with guided proposals" test-summary line.

### 7.5 Reproduce the actual numbers, both ways, fresh
This is the most important check in this whole document — it directly
underlies the session's final claim that the fixed kernel converges to a
better and more stable likelihood than upstream's.

**Against Aaron's repo** (patch only the two blocking bugs — missing `end`,
missing `Distributions` — nothing else):
```
cd ../PhyloPOMP.jl_aarons
julia --project=test -e 'using Pkg; Pkg.add("Distributions")'  # if not already done
cd .. && julia --project=PhyloPOMP.jl_aarons/test -e '
using Test, Crayons
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
src = read("PhyloPOMP.jl_aarons/test/mers_guided.jl", String)
Base.include_string(Main, src * "\nend\n", "mers_guided.jl")'
```
This takes ~5 minutes (170 mif iterations, no test assertions in the
original file — it will report "Total 0", which is expected and not a bug in
your reproduction). Instrument it (add `println("... = ", logLik(pf))` etc.
after each `@time` line) to actually see the numbers — do not accept a
silent completion as evidence of anything.

**Against our repo:**
```
cd PhyloPOMP.jl/test && RUN_HEAVY_TESTS=yes julia --project=. -e '
using Test, Crayons
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
include("mers_guided.jl")'
```
Instrument similarly.

Expected pattern (session's own numbers, for comparison — do not treat these
as ground truth, re-derive your own):

| stage | ours (fixed) | Aaron's (unfixed, patched only to run) |
|---|---|---|
| initial pfilter | ≈ −Inf or very poor | poor but finite |
| after mif, 120 iter | markedly better | better but worse than ours |
| after mif, 170 iter | best of the three checkpoints, ours | worse than its own 120-iter checkpoint |

CONCERN if your reproduction disagrees materially with this pattern — `mif`
is stochastic (LogNormal perturbations, no seed control after the first
`pfilter`), so exact numbers will differ, but the *qualitative* pattern
(ours ending substantially less negative) should hold across a few reruns.
FAIL if ours is *worse* than Aaron's unfixed version in a majority of reruns.

## 8. Task 7 — Compiler-layer claims (lower priority, verify if time allows)

```
cd PhyloPOMP.jl/test && RUN_HEAVY_TESTS=no julia --project=. -e '
using Test, Crayons
const h1 = crayon"bold blue"; const h2 = s -> crayon"!bold light_yellow"("- "*s)
include("kli_seir_compiled_test.jl"); include("kli_mers_compiled_test.jl")'
```
PASS iff both Gate-5 testsets report bit-exact agreement (worst `|Δll|` on
the order of `1e-14`, not merely "close"). This tests a claim that predates
this session (M08/M09) and should already be green; a regression here would
mean the POMP 0.10 bump broke something the handoff didn't check.

## 9. Task 8 — Git hygiene and scope discipline

```
git status --short
```

### 9.1 Nothing committed
PASS iff `git log --oneline -3` shows no new commits from this session (the
handoff explicitly says nothing was committed — verify, don't assume).

### 9.2 Modified-file list matches the handoff exactly
Cross-check every file in `git status --short`'s `M` list against
`handoffs/M14_upstream_guided_sync.md`'s "Files touched" section. FAIL if any
modified file is unexplained.

### 9.3 Untracked files — scope check
As of this review, `git status --short` shows several untracked files.
Classify each:
- `handoffs/M14_upstream_guided_sync.md`, `mers_guided_part1.md` (note
  prepended) — IN SCOPE, expected.
- `docs/compiler/guided_kernel_assessment.{tex,pdf}` — IN SCOPE (addendum
  added to a pre-existing untracked draft; confirm via git blame / session
  transcript that this file predates the session, and this session only
  edited its content, not created it from nothing).
- `docs/compiler/phylopomp_talk.*` — check timestamps; if these predate this
  session's start, they are NOT part of this work — flag but do not review
  their content as if they were.
- `docs/pipeline_plan.{tex,pdf}` — **explicitly flagged during the session as
  NOT created by this work** (timestamps didn't match any agent's activity).
  Confirm this file is still untouched by anything claiming to be part of
  M14; do not let it get swept into the "reviewed and approved" set.
- `docs/session_report_2026-09-19.{tex,pdf}`, `docs/tutorial/` — **appeared
  after the M14 work was reported complete, from a source this session did
  not identify.** Treat as explicitly OUT OF SCOPE for this review. Do not
  review, approve, or modify their content under this protocol. If asked to
  review them, that is a different, separate task — say so rather than
  silently including them.

### 9.4 No stray temp files
```
find /tmp/claude-1000 -newer handoffs/M14_upstream_guided_sync.md -name '*.jl' 2>/dev/null | grep -v tasks
```
This should be empty or only contain files this review protocol itself
created. Scratch scripts from the session's own investigation were not
preserved in the repo (by design) — verify none leaked into a tracked
directory:
```
git status --short | grep -v '^ M\|^??' # should be empty; catches staged-but-uncommitted anomalies
find src test docs handoffs -newer .git/HEAD -name '*.jl.orig' -o -newer .git/HEAD -name '*~' 2>/dev/null
```

## 10. Final report

For each of the 9 tasks above, give one line: PASS / FAIL / CONCERN, plus the
single most important piece of evidence. Then, separately, list every
CONCERN and FAIL with enough detail (command + actual output) that someone
who has not read this protocol can act on it without re-running everything
themselves. Do not average a FAIL against several PASSes into an overall
"looks good" — a single unverified likelihood-affecting bug fix (Task 3)
outweighs ten clean documentation diffs.
