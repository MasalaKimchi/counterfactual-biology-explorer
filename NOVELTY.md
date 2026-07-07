# Novelty — what is actually new, and how to make it newer

*Companion to `README.md` (framing) and `ASSESSMENT.md` (technical review). This file
answers three questions an expert judge will ask: **what is the novelty, stated as a
precise delta against the right prior art; how do we sharpen it; and where does it have
real-world impact.** The current docs position only against Mogrify and CellOracle —
that is the wrong comparison class and it undersells the contribution. The sharper
positioning is below.*

---

## 0. The one sentence (say this, not the paragraph in CLAUDE.md §0)

> **We turn a genome-scale CRISPRi screen into a feasibility oracle for cell-state
> engineering: given a desired transcriptional shift, we return a *verdict* — reachable
> by some combination of knockdowns, or provably not — together with the smallest
> knockdown set that reaches it and a certificate for what it cannot.**

The word doing the work is **verdict**. Everything else in the field returns a *ranking*
(which gene is most associated) or a *prediction* (what will this perturbation do). We
return a **decision with a proof attached**: reachable, or here is the direction no
knockdown can produce.

![Reachability as convex geometry: a feasibility verdict, not a ranking]({{artifact:art_152ac7f5-f20b-4521-a951-b98e6e0ec276}})

*Left: the target shift `d` lies inside the cone spanned by the (non-negative) knockdown
effect vectors — reachable, and the mixing weights are the minimal recipe. Right: `d`
lies outside the cone — no non-negative combination reaches it; the residual is a
direction that requires gene **activation**, and the separating hyperplane is the
infeasibility certificate. That right-hand panel is the output no existing method
produces.*

---

## 1. Why the current framing "does not make much sense" — and the fix

The paragraph in `CLAUDE.md §0` fails for a specific, fixable reason: **it names the
machinery (convex cone, non-negative combination, residual) before it establishes the
question the machinery answers.** A reader hits "convex cone" with no anchor for *why a
cone*, so the geometry reads as decoration on top of a ranking. It also buries the single
most novel output — the *provably-outside* verdict — inside a subordinate clause.

Rewrite it as **question → why the geometry is forced → what falls out**:

> A genome-scale CRISPRi screen measures, for every gene, the causal effect of knocking
> it down on the whole transcriptome. Read each of those measured effects as a vector,
> and a natural engineering question appears: **to move a cell from state A toward a
> target state B, which knockdowns — and in what combination — get closest, and is it
> even possible?**
>
> Two facts about CRISPRi make this question geometric. Knockdown only *removes* function
> (loss-of-function), and you can apply several knockdowns together but you cannot apply a
> "negative knockdown." So the reachable set of transcriptome shifts is exactly the
> **non-negative combinations of the measured effect vectors** — a convex cone. The
> target either lies inside that cone (**reachable** — and the combination weights are the
> recipe) or outside it (**provably unreachable by knockdown alone** — and the leftover
> direction tells you what would have to be *activated* instead).
>
> That verdict is the product. It is not a differential-expression list (what differs
> between A and B) and not a similarity score (which single gene looks most like the
> target). It is a **reachability decision with a residual and a confidence score**, and
> either answer — reachable or not — is a result you can take to the bench.

Notice this version leads with the *question* and earns "convex cone" by first stating
the two facts that force it, rather than opening on the machinery. ("Dictionary" is fine
as a technical term — it is standard sparse-coding vocabulary for a matrix of atoms — but
it belongs in the method section, not the opening sentence a non-specialist reads first.)

---

## 2. The novelty, stated as a precise delta against the RIGHT prior art

Novelty is a claim about a neighbourhood. The repo currently positions against **Mogrify**
and **CellOracle** — both about *reprogramming factor discovery* — and stops there. That
misses the field your method actually lives in, and against which the delta is cleanest:
the **network-control theory of cell fate**, a decade-old, active line of work.

### 2a. The comparison class you are missing: control of cell fate

This field asks *your exact question* — what is the minimal set of interventions that
drives a cell to a target fate — and has a mature toolbox for it:

| Prior method | What it does | What it controls | The basis it needs |
|---|---|---|---|
| **Zañudo & Albert 2015** (*PLoS Comput Biol*) — stable-motif control | Minimal node interventions that switch a Boolean network's attractor | discrete attractors | a hand-built **Boolean model** |
| **NETISCE 2022** (*npj Syst Biol Appl*) — feedback-vertex-set control | Ranks intervention sets that steer a signalling network to a desired steady state | steady states | an **inferred/curated network** + simulation |
| **CellOracle 2023** (*Nature*) | Simulates KO/overexpression through an inferred GRN into a "vector map of transitions" | cell-identity transitions | a **GRN inferred from WT scRNA-seq + ATAC** |
| **"Reverse control of cell-fate" 2024** (*npj Syst Biol Appl*) | Derives control targets from cell-fate landscape structure | fate attractors | a **model of the fate landscape** |
| **Mogrify 2016** (*Nat Genet*) | Predicts TF set for a transdifferentiation | cell-type conversion | expression + **regulatory-network** heuristic |

**Every one of these runs on an *inferred or hand-built model* of the regulatory
network.** That is the field's structural bottleneck, and it is stated openly in the
literature: the control set you get out is only as trustworthy as the network you put in,
and GRN inference from observational data is famously underdetermined. These methods
answer "what should we perturb" by *simulating a model* of what perturbations would do.

### 2b. Your delta — three things, in priority order

**(1) The effect vectors are measured, not modelled — this collapses the field's
bottleneck.** Your dictionary rows are the *empirically observed* genome-scale CRISPRi
effects on the real transcriptome. You are not inferring what a knockdown does and then
controlling the inference; you are composing measurements. This is the difference between
"control the model and hope it transfers" and "compose interventions whose effects were
directly observed in the target cell type." *This is the headline — lead with it.*

**(2) A convex-cone reachability *verdict with a certificate*, not a ranked set.** The
control-theory methods return ranked/heuristic intervention sets. None of them return a
*provable infeasibility result*. Because your reachable set is literally a convex cone,
"unreachable" is not "we didn't find a good set" — it is a theorem: the target has a
component in the orthogonal complement of the non-negative span, and the separating
hyperplane is a checkable certificate. **"Provably outside the knockdown cone → this
transition requires activation" is an output no method in the table above produces.** A
negative result becomes a positive, falsifiable claim.

**(3) Held-out-gene validity + a screen-native confidence decomposition.** You validate
by fitting weights on a random half of the target signature's genes and scoring alignment
on the held-out half — a generalisation test that cannot be inflated by overfitting, and
that the association/ranking framings structurally cannot run. Confidence decomposes into
the screen's own reproducibility metrics (cross-guide, cross-donor), stability-selection
frequency, and orthogonal literature/genetics evidence.

### 2c. What you must NOT claim as novel (credibility depends on this)

- **Not the Th1/Th2 or aging regulators.** Zhu et al. already report these. Recovering
  GATA3 is a *positive control*, not a discovery. (`ASSESSMENT.md §1` is already honest
  about this — keep it that way.)
- **Not "combining perturbation effects linearly."** Additive composition of effect
  vectors is standard; the 2025 *Nat Methods* baseline paper (Ahlmann-Eltze et al.) is
  your friend here — it licenses the linear model as the right primary model, but linear
  composition itself is not the novel part.
- **Not "convex optimisation on omics."** NNLS/LASSO on expression data is old. The novel
  composite is *measured LOF effect vectors + non-negativity-as-biology + a feasibility
  verdict with a certificate*, applied to cell-state engineering.

### 2d. The defensible novelty sentence (drop-in)

> Prior cell-fate control methods (stable-motif control, NETISCE, CellOracle, reverse
> control) compute intervention sets by simulating an **inferred model** of the regulatory
> network. We instead compose the **measured** genome-scale CRISPRi effect vectors
> directly, and — because knockdown is loss-of-function and combinations are non-negative
> — cast reachability as membership in a convex cone. This yields what model-simulation
> approaches do not: a **provable reachable / outside-the-cone verdict with an
> infeasibility certificate**, validated by held-out-gene generalisation. We claim novelty
> for this measured-effect, certificate-carrying reachability decision — not for the
> regulator lists, which the source screen already reports.

---

## 3. How to make the novelty *stronger* (ranked by value-per-effort)

These are ways to widen the delta, not just defend it. ★ = do this for the hackathon;
the rest are the "future work" slide.

**★ (a) Ship the infeasibility certificate as a first-class output, not a residual
number.** Right now "outside the cone" is planned as a scalar residual. Make it
*constructive*: when the target is unreachable, return the **specific genes in the
residual direction that would need activation** — i.e. the dual variables / the
separating hyperplane's heavy coordinates. That converts "not reachable" into "not
reachable by knockdown; here are the 3 genes a CRISPRa arm would need to test," which is
a concrete, wet-lab-actionable hypothesis. **This is the single highest-value addition**
— it is the output that most sharply distinguishes you from every ranking method, and it
directly motivates the Schmidt 2022 CRISPRa validation you already cite.

**★ (b) Report the reachability *spectrum* as the headline figure, with both nulls.**
Alignment-vs-k against a shuffled-target null *and* a random-perturbation null (already in
the roadmap). This is what makes the verdict honest — the null band is the
achievable-by-chance floor at every sparsity level — and it is the figure that survives
expert scrutiny. Without it, any single cosine is dismissable as overfitting.

**★ (c) Make non-negativity earn its keep with a directional decomposition.** Split the
target `d` into its *in-cone* (knockdown-reachable) and *out-of-cone* (activation-
requiring) components and report the fraction of the target norm in each. "72% of the
Th2→Th1 shift is reachable by knockdown; 28% requires activation, concentrated in
{TBX21, STAT4, ...}" is a far richer, more falsifiable claim than a single verdict, and
it is a direct consequence of the geometry you already have.

**(d) Turn the additivity caveat into a model feature, then a designed experiment.**
Penalise co-selecting genes with collinear effect vectors (shared downstream program →
predicted sub-additivity), then rank the single most informative *next* experiment as the
pair whose additive-vs-epistatic prediction is most decision-relevant. This reframes
"combinations are untested" from a limitation into the product: a prioritised experiment
queue. (Already sketched in `ROADMAP.md` item 21 — promote it.)

**(e) Condition-resolved cones.** You have Rest / Stim8hr / Stim48hr. Ask whether a target
is reachable in one condition but not another — "reachable only in activated T cells" is a
biologically meaningful and novel kind of statement, and it uses data you already have.

**(f) Cross-cone transfer as the generality claim.** Nothing in the method is
T-cell-specific. Demonstrating the *same* machinery on a second public Perturb-seq atlas
(even a small one) upgrades the contribution from "a T-cell result" to "a general operator
for perturbation atlases." This is the strongest *scientific* claim and the natural
future-work headline.

**(g) Bound the cone's resolution honestly.** With ~2,000 perturbations spanning a
10,000-gene space, the cone is low-dimensional relative to the ambient space, so *most*
random targets will look "partially outside." State this — it is why the shuffled-target
null (b) is essential, and acknowledging it pre-empts the sharpest reviewer objection.

---

## 4. Real-world impact — where a feasibility oracle is actually useful

The impact is not "we found a drug target." It is that **a reachable/unreachable verdict
changes what a lab does next**, and it does so before any wet-lab spend. Concrete uses:

### 4a. Cell-therapy and immuno-engineering target triage
The featured dataset is CD4+ T cells; the immediate application is **engineering T-cell
state for therapy** (CAR-T persistence, Treg stability, Th-bias correction). A programme
deciding which knockdowns to put into an arrayed or in-vivo screen faces a combinatorial
budget. A tool that says *"this state shift is reachable by this 3-gene knockdown set (do
this), and this other shift is NOT reachable by knockdown — don't waste a CRISPRi arm on
it, it needs activation"* is **directly decision-relevant** and saves the most expensive
resource in the pipeline: screen slots and animal cohorts.

### 4b. Autoimmune / allergic disease — the Gladstone-prize hook
The local data already links perturbation clusters to **17 immune-mediated diseases with
185 significant GWAS-gene enrichments** (Crohn's, IBD, asthma, atopic eczema, psoriasis,
MS, T1D, lupus, RA, ...). Th1/Th2 balance is a validated therapeutic axis — Th2-skew
drives allergy/atopy, Th1/Th17 drives autoimmunity. A method that nominates the *minimal,
ranked, confidence-scored* knockdowns to rebalance that axis — **and states when a
rebalancing is not achievable by knockdown at all** — is a target-selection instrument for
immunotherapy, with the disease link already in the data rather than asserted.

### 4c. The "don't run this experiment" value — the underrated one
Most target-nomination tools only ever say *go*. The distinctive value of a feasibility
oracle is that it also says **stop**: a provable "outside the knockdown cone" verdict tells
a team that a loss-of-function screen *cannot* reach their goal and that they need an
activation modality (CRISPRa, an agonist) instead. **Negative results that redirect
resources are where a decision tool earns its cost** — and this is the use case no ranking
or association method can serve.

### 4d. General instrument for the Perturb-seq era
Perturb-seq atlases are proliferating. Any atlas + any target state (disease→healthy,
aged→young, cell-type A→B) plugs into the same cone machinery. The larger claim is a
**reusable reachability operator for perturbation atlases** — the transferability, not the
one T-cell result, is the scientific contribution.

### 4e. Honesty guardrail (state all three, plainly)
Outputs are **hypotheses for wet-lab testing**, not validated targets. CRISPRi is
loss-of-function only. Multi-gene sets assume additivity (no epistasis), and matching a
transcriptional signature is not proof of functional rescue. The honest framing is what
makes the impact claim credible — an oracle that over-promises is worthless to a lab.

---

## 5. One-line summary for the pitch

> Everyone else ranks genes or predicts what a perturbation will do. **We answer a
> yes/no engineering question with a proof: can you get this cell to that state by
> knockdowns, what's the smallest set, and — when you can't — what must you activate
> instead.** Measured effects, convex geometry, a certificate either way.
