# Review brief: Milestone 4 — the BIRTH branch had the same bug class as the SAMPLE fix

**Audience:** a reviewing agent that has not seen the implementation session.
Assumes `check_milestone1.md`–`check_milestone3.md` read first.

**Branch:** `atpabuser-devel`. **Status:** bug found, fixed, and independently
re-validated against R. Not yet committed at the time of writing.

## 1. How this was found

The user asked for a rigorous, evidence-based case that the simulator is
correct, explicitly requesting deeper ("Opus-level") review of the existing
work rather than a restatement of it. That review re-read `src/simulate.jl`
against R's actual C++ kernel (`~/Desktop/phylopomp-fork/src/master.h`) —
not just re-running the existing tests — and found that `apply_event!`'s
`BIRTH` branch has **the same bug class** Milestone 3 found and fixed for
`SAMPLE`: the node representing the divergence point is stamped with the
*previous* event's time on the parent lineage, not the birth's own firing
time.

This was independently re-derived and confirmed in-session (not taken on
the reviewing agent's word alone): a direct instrumentation of `simulate`
showed 96.5% of internal two-child nodes strictly precede both of their
children in time by a nonzero, model-dependent lag (mean ≈1.2 time units
at the test parameters) — the exact signature described.

## 2. Why it's real, and why nothing in Milestones 1–3 caught it

**The prior code** (`apply_event!`, `BIRTH` branch):
```julia
remove!(inv,d,b)
c1 = push_node!(G,t,Node,b)      # child 1, correctly stamped t
push!(cur.children,c1)
add!(inv,d,c1)
for j ∈ ev.into
    cj = push_node!(G,t,Node,b)  # child 2, correctly stamped t
    push!(cur.children,cj)
    add!(inv,j,cj)
end
```
Both new children ARE correctly stamped at the firing time `t` — but they
are attached as children of `cur`, whose own `slate` is the time of the
lineage's *previous* event (`t_prev < t`). Since `newick.jl` writes branch
length as `child.slate - parent.slate`, and both the Julia filters and R
read an internal node's own `slate` as its divergence time, the fork is
recorded at `t_prev`, not `t` — the split is displaced backwards in time by
the elapsed time since `cur`'s own creation.

**Why every existing check missed it:**
- The self-consistency invariant (`length(inv[d]) == x[demes[d]]`) is a
  *count*, blind to which node carries which timestamp.
- The Newick/CBLV round-trip checks node/sample *counts* and (in
  `seir_simulate.jl`) `times(g)` as an unordered vector compared against
  itself post-round-trip — self-consistent even when every time is wrong,
  because both sides of that comparison inherit the same error.
- The filter-agreement headline test (`NaiveSEIR.filter_pomp` /
  `SoftMERS.filter_pomp` returning finite `logLik`) checks structural
  legality (node degrees) and rate-compatibility of the branching pattern,
  not *when* a branch point falls. A displaced-but-still-bifurcating,
  still rate-compatible tree still parses and still filters.
- The Milestone-3 cross-validation's three statistics (extinction fraction,
  `nsample`, first-sample-time) are a count, a count, and the timestamp of
  a `SAMPLE` node specifically — `SAMPLE` was exactly the branch that had
  already been fixed. None of them touch an internal `Node`'s own time.

This is the direct, sharper version of the reviewer-focus note already
logged in `check_milestone3.md` §7.1 ("re-derive independently that the
interposed-`Node` handling is correct, rather than trusting this
document's trace") — re-deriving it one step further is exactly what
surfaced this.

**Measured impact** (this session's own instrumentation, independent of the
subagent's numbers): 96.5% of internal branch points affected, mean lag
1.21 time units, max 7.93, at the same parameters used throughout this
milestone series. Every zero-length edge this produced also broke the
Newick round trip about as often (any lineage whose *next* event happened
to be a birth needed a same-instant continuation node, which the writer
happily emits and `parse_newick`'s `clip_zlb!` — correctly, by design —
refuses to collapse cleanly under a `Sample` parent, and merges children
into a polytomy under a plain-`Node` parent when the two branches beneath
it hadn't both survived pruning as exactly one child each).

## 3. The fix

Read R's `genealogy_t::birth` (`src/master.h`) directly:
```cpp
ball_t* birth (ball_t* a, slate_t t, name_t d) {
  time() = t;
  node_t *p = make_node(a->deme());   // ONE fresh node, stamped t
  ball_t *b = new ball_t(p,p->uniq,black,d);
  p->insert(b);                       // the new daughter lineage
  p->slate = time();
  add(p,a);                           // AND the continuing lineage
  return b;                           // -- both now open at p
}
```
One new node at the firing time, shared by **both** resulting lineages;
neither gets a node of its own until *its own* next event fires. Ported
directly (`src/simulate.jl`, `apply_event!`):
```julia
if ev.type == BIRTH
    remove!(inv,d,b)
    p = push_node!(G,t,Node,b)
    push!(cur.children,p)
    add!(inv,d,p)
    for j ∈ ev.into
        add!(inv,j,p)
    end
```
Consequence worth stating explicitly: the same node name `p` can now
appear in **two different demes'** inventory vectors simultaneously
(once for the continuing lineage, once per new lineage) — this is
intentional and matches R's "one node, two balls" representation exactly.
`remove!`/`add!` are deme-scoped, so this doesn't create any ambiguity:
each deme's vector independently tracks its own membership, and `p` simply
represents "two logical lineages currently colocated at one physical tree
position," resolved into two separate children only once each lineage has
its own next event. The self-consistency invariant continues to hold by
the same count-preserving argument as before (verified: full suite green).

A pleasant side effect: this also makes the `Sample`-typed-`cur` special
case from Milestone 3 (interposing an extra plain `Node` so a `Sample`
node doesn't directly gain two children) unnecessary — a `Sample`-typed
`cur` now gets exactly one new child (`p`) either way, same as any other
`cur`, so the special case was removed. The code got shorter as a result of
being made more correct, which is a good sign, not a coincidence: the two
bugs were the same shape, and fixing the general case subsumed the
band-aid over the specific one.

## 4. Validation

**Direct, R-independent**: 400 non-extinct SEIR draws, same parameters as
throughout — zero zero-length edges (was up to 1018 across 289/312 trees),
zero Newick round-trip node-count mismatches (was 276/312), max node
degree 2 everywhere (`Root` 1, internal `Node` 2, `Sample` ≤1, confirmed
directly, not merely absence-of-assertion-failure).

**Cross-validation against R, before vs. after** (800 draws/side, same
seeded methodology as `scripts/seir_crossvalidate.{R,jl}`, extended with
two new statistics this fix specifically targets):

| statistic | before fix (KS D, p) | after fix (KS D, p) |
|---|---|---|
| `nsample` | 0.036, 0.54 (already fine) | 0.036, 0.67 |
| mean internal-node time | **not measured by the committed script; independently confirmed broken (p≈0) by the reviewing subagent** | **0.031, 0.94** |
| total branch length | **same — inflated ~15% (356 vs 308) per subagent measurement** | **0.049, 0.49 (means: 309.12 vs 309.13)** |
| number of internal nodes | 0.029, 0.96 | 0.029, 0.96 |

The two statistics the bug actually broke go from a subagent-reported
`p ≈ 0` / ~15%-inflated mean to `p > 0.4` and means agreeing to three
significant figures, while `nsample` — unaffected by this bug either way —
stays exactly where it was. That pairing (the broken statistics move, the
unrelated one doesn't) is the direct evidence the fix addresses the actual
defect rather than perturbing the process some other way.

**Full regression suite**: `RUN_HEAVY_TESTS=no julia --project=.
test/runtests.jl` → **21266/21266 passing** (this count includes a large,
unrelated "KLI compiled filter" test track added to the repo by other work
during this session; the simulator's own acceptance tests,
`test/{seir,mers}_simulate.jl`, are 41/41 as before — unchanged in count
because none of them assert on internal-node timing, which is exactly the
gap `check_milestone3.md` §7.2's suggested follow-up would have closed).

## 5. What is now believed correct, and what still isn't

**Now correct, confirmed by direct measurement, not merely by the absence
of a failing assertion:** internal branch-point times, total branch
length, coalescent-interval structure, ranked topology — the statistics
that were wrong are the statistics that were checked, and they now agree
with R to within sampling noise.

**Still open**, unchanged from `check_milestone3.md` §7 and now sharpened
by this round:
- The class of bug found twice in one model (SEIR, the only one with any
  cross-language validation) has never been checked for MERS at all —
  MERS's sampling is always destructive, so it never exercised the
  `Sample`-then-`BIRTH` path either bug lived in, but MERS's four
  transmission events are still `BIRTH`s, subject to the same node-sharing
  logic, and have no independent time-based validation.
- The concrete, cheap fix this session did *not* add: a permanent
  unit-level regression asserting every emitted node's `slate` equals the
  firing time of the event that created it, for every event type, not
  just retrospectively via distributional comparison. `check_milestone3.md`
  §7.2 asked for this after the first bug; it would have caught this one
  at Milestone-1 speed instead of requiring a second, deeper review pass.
  Still not done — flagged again, more urgently, here.
- Multi-root/forest and piped continuation remain untested/unimplemented,
  as before.

## 6. Process note

This bug was found because the review was asked to go deeper than
confirming the existing test suite still passes — it re-read primary
source (R's actual C++, not a description of it) and re-derived the
correctness argument from scratch rather than accepting the prior
milestones' conclusions. The prior milestones' own conclusions were not
wrong about what they checked; they were incomplete about what they didn't
check, and said so explicitly (`check_milestone3.md` §7.1 named the exact
spot). Worth internalizing for any future round: "the tests pass" and "the
distributional comparison didn't reject" are claims about the *statistics
tested*, not about correctness in general, and the gap between them is
exactly where bugs like this hide.
