# Real-world impact — why a reachability oracle is useful, and to whom

*Companion to `NOVELTY.md`. Where NOVELTY explains what is scientifically new, this file
makes the case an industry scientist would need to hear: the problem is large and
expensive, the specific failure mode this method attacks is the dominant one, and the
"measured, not inferred" design is aligned with the only two things shown to move
drug-development odds. Numbers are cited inline; see the provenance note at the end for
which were read from the primary source this session versus confirmed at citation level.*

---

## 1. The problem is failure, and failure is the cost

Drug development is not primarily expensive because success is expensive — it is expensive
because **failure is the norm and it is paid for late.** The canonical figures:

- **~90% of clinical programs never reach approval.** Only about 10% of clinical
  programmes eventually receive approval (Minikel et al., *Nature* 2024; the BIO/Thomas
  2006–2015 analysis puts the cumulative rate at 9.6%). The cost of drug discovery and
  development is *driven primarily by failure*.
- **The failure concentrates at Phase II, and its cause is efficacy.** Phase II has the
  lowest transition rate of any phase — **30.7%**, versus 63.2% (Phase I) and 58.1%
  (Phase III) in the BIO dataset (Thomas et al. 2016). And **lack of efficacy is the
  reason for approximately half of all Phase II and Phase III failures** (Harrison, *Nat
  Rev Drug Discov* 2016, analysing 174 attrition events with a stated reason over
  2013–2015) — i.e. wrong-target / wrong-mechanism, not safety or operations, is the
  dominant late-stage failure mode.
- **Autoimmune indications — this dataset's home turf — sit near the bottom.** Wong, Siah
  & Lo (*Biostatistics* 2019, 406,038 trial records) put autoimmune candidates at roughly
  **15% probability of success**, well below the cross-indication average.

**The implication for a computational method:** the highest-leverage place to intervene is
*target selection* — the decision, made years before Phase II, of *which* gene to go
after. A wrong target is an efficacy failure waiting to happen, and it is only discovered
after the most expensive experiments have run.

---

## 2. Two things move the odds — and both are "measured, not inferred"

The literature is unusually clear about what *reduces* efficacy failure, and both levers
point the same way as this project's core design choice.

### 2a. Human genetic support ~2.6× the odds of approval
Drug mechanisms with human genetic support are **2.6 times more likely** to succeed from
Phase I to approval than those without (Minikel et al., *Nature* 2024); the earlier
estimate (Nelson et al., *Nat Genet* 2015) was ~2×. Critically, the authors frame genetic
support not as a gate but as **an enrichment signal and a probabilistic tool for portfolio
prioritisation.** That is exactly the register a reachability oracle operates in: it does
not promise a target works — it *re-weights the portfolio* toward mechanisms with
independent, measured support.

*This project's disease layer already carries this signal:* the local autoimmune-
enrichment table links perturbation clusters to **17 immune-mediated diseases with 185
significant GWAS-gene enrichments (FDR < 0.05)**. A nomination that is *both*
knockdown-reachable *and* genetically supported for the target disease is precisely the
"high relative-success" quadrant Minikel describes.

### 2b. Inferred preclinical claims don't reproduce — measured effects are the antidote
The reason "wrong target" is so common upstream is that the published preclinical
literature it draws on is not reliable. The two landmark industry audits:

- **Amgen: 6 of 53 (11%)** landmark preclinical cancer papers could be reproduced (Begley
  & Ellis, *Nature* 2012).
- **Bayer: ~25%** of target-validation projects reproduced cleanly (Prinz et al., *Nat Rev
  Drug Discov* 2011).

This is the crux of the "measured, not inferred" argument in `NOVELTY.md §2b`, stated in
business terms: **the field's target ideas are largely built on inferred networks and
individually-published effects that fail to replicate.** A method whose dictionary rows
are *directly measured genome-scale CRISPRi effects in the relevant primary human cell
type* is drawing on exactly the kind of systematic, internally-reproducible measurement
(cross-guide, cross-donor) that the inferred literature lacks. That is not a
nice-to-have — it is the specific defect that makes 90% of programs fail.

**The one-sentence pharma pitch:** *the two interventions shown to cut efficacy failure —
human-genetic support and measured (not inferred) target effects — are the two things this
method is built on. It scores knockdown nominations by measured effect and cross-references
them to disease genetics, in the primary cell type where the biology actually happens.*

---

## 3. Where it is useful — concrete users and decisions

### 3a. Target triage for T-cell engineering & cell therapy (the nearest application)
The featured dataset is primary human CD4+ T cells, and **T-cell state engineering is a
market in the middle of a step-change.** CAR-T, proven in oncology, is moving into
autoimmunity fast: as of December 2024 there were **116 clinical trials** evaluating CAR-T
against autoimmune conditions (Clinical & Translational Science, 2025), and trials in
autoimmune rheumatic disease alone **peaked at 25 in 2024** (Frontiers in Immunology,
2025), with early SLE cohorts reaching durable, steroid-free remission.

A programme deciding *which knockdowns* to engineer into a T-cell product — for
persistence, for Treg stability, for Th-bias correction — faces a combinatorial design
space and a fixed screening budget. **This is the decision the reachability oracle is
built for:** "this state shift is reachable by this ranked 3-gene knockdown set — put it in
the arrayed screen" and, just as valuably, "this shift is *not* reachable by knockdown —
don't spend a CRISPRi arm on it." The users are cell-therapy target-selection teams and
academic immuno-engineering labs; the decision is which constructs enter a screen that
costs weeks and animal cohorts per arm.

### 3b. The "don't run this experiment" verdict — the underrated value
Every ranking or association tool only ever says *go*. The distinctive output here is a
**provable stop**: when the target lies outside the knockdown cone, the tool returns not
just "no" but the **specific genes that would have to be *activated* instead** (the
constructive certificate in `reachability.py`). For a lab, that converts a dead end into a
redirected experiment: *stop the loss-of-function screen, switch to CRISPRa / an agonist,
and here are the 3 genes to test first.* Given that late efficacy failure is the dominant
cost (§1), **a credible early "this modality cannot reach the goal" is worth more than
another optimistic ranking** — it moves resources before they are sunk.

### 3c. Portfolio prioritisation across a Perturb-seq atlas
Framed the way Minikel frames genetic support — a probabilistic prioritisation signal, not
a gate — the reachability verdict + confidence decomposition is a **portfolio-ranking
instrument.** For any (target state, disease) pair, it returns a reachable/outside verdict,
a minimal ranked set, a screen-native confidence score, and a disease-genetics
cross-reference. A discovery team can rank a whole slate of desired cell-state changes by
*how reachable they are and how well-supported the drivers are* — before committing a
single wet-lab FTE.

### 3d. A reusable operator for the Perturb-seq era (the general claim)
Nothing in the method is T-cell-specific. Perturb-seq atlases are proliferating across
tissues and disease models; any atlas + any target state (disease→healthy, aged→young,
cell-type A→B) plugs into the same convex-cone machinery. The durable contribution is a
**general reachability operator for perturbation atlases** — the transferability, not the
one T-cell result, is the scientific and commercial headline.

---

## 4. Quantified value proposition (order-of-magnitude, honestly framed)

The value is not a claimed hit rate — it is **shifting where a program spends and where it
stops.** In the terms of §1–§2:

| Lever this method pulls | The number it speaks to | Direction |
|---|---|---|
| Prioritise measured over inferred target effects | 11–25% preclinical reproducibility (Amgen/Bayer) | attacks the root cause of wrong-target efficacy failure |
| Cross-reference nominations to disease genetics | 2.6× approval odds with genetic support (Minikel 2024) | moves the portfolio toward the high-success quadrant |
| Provable "unreachable by knockdown" verdicts | ~48% of Phase II failures are efficacy; autoimmune PoS ~15% | redirects modality *before* the expensive phase |
| Rank minimal knockdown sets for a screen | screen slots / animal cohorts per arm | allocates the scarcest preclinical resource |

**The framing that survives scrutiny:** this is a **hypothesis-prioritisation and
experiment-triage instrument**, not a target-validation engine. It tells a team where to
look first and where not to look at all, using measured effects and disease genetics — the
two signals with published evidence of improving odds. It does not replace the wet lab; it
decides what the wet lab does next.

---

## 5. Honesty guardrails (state these plainly — they are what make the impact credible)

- **Nominations are hypotheses for wet-lab testing, not validated targets.** A
  transcriptomic-signature match is not proof of functional rescue, let alone clinical
  efficacy.
- **CRISPRi is loss-of-function only.** The tool can nominate knockdowns and *flag* what
  needs activation; it cannot test the activation hypotheses itself (that is the CRISPRa
  validation arm, e.g. Schmidt et al. 2022).
- **Multi-gene sets assume additivity (no epistasis).** Every combination is an
  extrapolation to be tested, and should be flagged as such.
- **One primary-cell system, four donors.** External validity across genotypes, tissues,
  and disease contexts is a claim to be earned, not assumed.

An oracle that over-promises is worthless to a drug-development team — the honest framing
above is precisely what makes the (real, large) upside believable to an expert.

---

## References

1. Minikel, Painter, Dong & Nelson. *Refining the impact of genetic evidence on clinical
   success.* **Nature** (2024). doi:10.1038/s41586-024-07316-0. — verified quotes: "only
   about 10% of clinical programmes eventually receiving approval"; "probability of
   success for drug mechanisms with genetic support is 2.6 times greater than those
   without."
2. Nelson et al. *The support of human genetic evidence for approved drug indications.*
   **Nat Genet** 47, 856–860 (2015). — original "genetic evidence doubles the success
   rate" (~2×) estimate, as restated in ref [1].
3. Wong, Siah & Lo. *Estimation of clinical trial success rates and related parameters.*
   **Biostatistics** 20(2):273–286 (2019). doi:10.1093/biostatistics/kxx069. — 406,038
   records over 21,143 compounds; autoimmune PoS ~15%; oncology 3.4%.
4. Thomas et al. *Clinical Development Success Rates 2006–2015.* BIO / Biomedtracker /
   Amplion (2016). — 9.6% cumulative success; Phase II transition 30.7% (vs Phase I 63.2%,
   Phase III 58.1%).
5. Harrison. *Phase II and phase III failures: 2013–2015.* **Nat Rev Drug Discov** 15,
   817–818 (2016). doi:10.1038/nrd.2016.184. — verified quote: "lack of efficacy is the
   reason for approximately half of all phase II and phase III failures" (174 failures
   with a stated reason).
6. Begley & Ellis. *Drug development: raise standards for preclinical cancer research.*
   **Nature** 483, 531–533 (2012). doi:10.1038/483531a. — Amgen reproduced 6 of 53 (11%)
   landmark papers.
7. Prinz, Schlange & Asadullah. *Believe it or not: how much can we rely on published data
   on potential drug targets?* **Nat Rev Drug Discov** 10, 712 (2011). — Bayer ~25%
   reproducibility.
8. CAR-T in autoimmunity: 116 trials as of Dec 2024 (Clinical & Translational Science,
   2025); autoimmune-rheumatic-disease trials peaked at 25 in 2024 (Frontiers in
   Immunology, 2025).

*Provenance note (honest about verification level). Read at sentence level from the
primary source this session: the 2.6× / ~10%-approval figures (ref 1), the ~50%-of-
failures-are-efficacy figure (ref 5, Harrison 2016), Amgen 6/53 (ref 6), and the CAR-T
counts of 116 and 25-in-2024 (ref 8). Consistent with well-known published figures but
confirmed here only at title/citation level, not re-read from the source text this
session: the 9.6% cumulative rate and the 30.7% / 63.2% / 58.1% phase-transition rates
(ref 4, BIO/Thomas 2016), the ~15% autoimmune PoS (ref 3), the ~2× Nelson 2015 estimate
(ref 2), and Bayer ~25% (ref 7). Confirm the latter set against the publisher record
before any external submission.*
