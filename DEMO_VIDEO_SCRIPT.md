# 3-Minute Demo Video — Shooting Script

**Project:** Counterfactual Biology Explorer / Cell-State Reachability
**Event:** Built with Claude — Life Sciences Hackathon (Research Track)
**Runtime target:** 2:55–3:00 · **Voiceover:** ~400 words @ ~150 wpm

---

## The logline (say this if you have one sentence)

> **Every AI model in biology predicts what a perturbation *will* do. We built the
> first one that tells you what's *impossible* — before you run the experiment.**

---

## The hook (recommended) — "the word no other model can say"

Open on the pain pharma feels in its bones (attrition), land the reframe
(*feasibility*, not prediction), and prove the novelty with the empty quadrant.

**Cold-open line (0:00):**
> "Nine of ten drug programs fail — most die in Phase Two, and the top cause is
> lack of efficacy. Today's AI predicts what a perturbation *will* do. We built
> the first tool that tells you what's *impossible* — before you run it."

### Three alternate hooks (swap the cold open if you prefer)

1. **The GPS line (metaphor-first, warmest):**
   *"A GPS doesn't just find the fastest route — it tells you when there's no road
   at all. We built that for cell engineering."*
   → Use as the mid-video explainer even if you don't open with it (Beat 3).
2. **The stop-loss line (finance-flavored, blunt):**
   *"Cell engineering has no stop-loss. You find out the target was unreachable
   after the trial fails. We move that verdict to before the screen."*
3. **The empty-quadrant line (novelty-first, for a methods-savvy panel):**
   *"We surveyed 91 methods across 15 years. Zero can tell you a cell state is out
   of reach. That quadrant was empty. Now it isn't."*

---

## Beat sheet / shot list

| # | Time | On screen (capture source) | Voiceover |
|---|------|----------------------------|-----------|
| 1 | 0:00–0:18 | **Attrition funnel** — screen-record `app/explorers/pharma_funnel.html`; let the "Provable 'unreachable' STOP" flag sit | "Nine of ten drug programs fail — most die in Phase Two, and the top cause is lack of efficacy. Today's AI predicts what a perturbation *will* do. We built the first tool that tells you what's *impossible* — before you run it." |
| 2 | 0:18–0:40 | Dataset beat — `fig_central_illustration.png` Panel A (CRISPRi screen → effect matrix), slow push-in | "Start with a real genome-scale CRISPRi screen — thirty-four thousand knockdowns, measured in primary human T cells. The usual question is *what's different* between a sick cell and a healthy one. We ask a harder one: is the healthy state even *reachable* — and if so, what's the smallest set of knockdowns that gets there?" |
| 3 | 0:40–1:05 | **The cone** — `app/explorers/reachability_explorer.html` Panel A (live), or Fig 1 Panel B | "Here's the idea. Each knockdown is a measured effect vector. Because knockdown only ever subtracts, the states you can reach form a *cone*. Reachability is one question: does your target live inside that cone? A GPS doesn't just find the fastest route — it tells you when there's no road at all. This is that, for cell state." |
| 4 | 1:05–1:33 | **The verdict + certificate** — `app/explorers/reachability_explorer.html` Panel B (signed decomposition), then Fig 1 Panel C activation certificate | "The answer isn't a similarity score — it's a falsifiable verdict. For the Th2-to-Th1 switch: *partly reachable*. Thirty-nine percent of the target is reachable by knockdown; twenty-five percent provably is not — it needs genes turned *up*. And the tool names the exact genes no knockdown mix can deliver. That's a wet-lab prediction you can go falsify tomorrow." |
| 5 | 1:33–2:00 | **Modality triage** — screen-record `app/explorers/pharma_triage.html`; hover IRF1, JAK2, ICOS | "Now the part pharma cares about. Take the hundred-and-two knockdowns the model says you actually need — and cross them against druggability and human genetics. Forty-four percent are hard-to-drug. Ten are clinical-grade today. IRF1 lights up: strong genetics, undruggable. JAK2 and ICOS: green-light. That's a target shortlist triaged by *feasibility* — not a ranked gene list." |
| 6 | 2:00–2:20 | **The differentiator** — `app/explorers/pharma_capability.html`; zoom the empty top-right quadrant, then the lone "this work" dot appears | "We surveyed ninety-one prior methods — fifteen years of them. Fourteen use measured data. *Zero* return a feasibility verdict. Zero can prove something is out of reach. This whole quadrant — measured *and* feasible — was empty. Now it isn't." |
| 7 | 2:20–2:42 | **Trust + transfer** — Fig 1 Panel D (K562 transfer), flash `app/explorers/causal_trust.html` | "And it travels. The same operator, unchanged, jumps from T cells to a completely different screen — different cell type, opposite direction — and still recovers a held-out target at cosine point-eight-eight. Every number ships with a null model, a confidence score, and a one-command reproduce." |
| 8 | 2:42–3:00 | **Close** — `app/index.html` narrative scroll (the six-chapter walkthrough that threads all seven explorers), end on title card + logline | "Seven interactive explorers, a full manuscript, a ninety-one-method survey — built in a week with Claude Science. We answer the question that saves the most money in drug discovery: before you run it — can you even get there?" |

---

## Clean voiceover script (for a narrator or TTS, ~400 words)

> Nine of ten drug programs fail — most die in Phase Two, and the top cause is lack
> of efficacy. Today's AI predicts what a perturbation *will* do. We built the first
> tool that tells you what's *impossible* — before you run it.
>
> Start with a real genome-scale CRISPRi screen — thirty-four thousand knockdowns,
> measured in primary human T cells. The usual question is *what's different* between
> a sick cell and a healthy one. We ask a harder one: is the healthy state even
> *reachable* — and if so, what's the smallest set of knockdowns that gets there?
>
> Here's the idea. Each knockdown is a measured effect vector. Because knockdown only
> ever subtracts, the states you can reach form a *cone*. Reachability is one question:
> does your target live inside that cone? A GPS doesn't just find the fastest route —
> it tells you when there's no road at all. This is that, for cell state.
>
> The answer isn't a similarity score — it's a falsifiable verdict. For the
> Th2-to-Th1 switch: *partly reachable*. Thirty-nine percent of the target is reachable
> by knockdown; twenty-five percent provably is not — it needs genes turned *up*. And
> the tool names the exact genes no knockdown mix can deliver. That's a wet-lab
> prediction you can go falsify tomorrow.
>
> Now the part pharma cares about. Take the hundred-and-two knockdowns the model says
> you actually need — and cross them against druggability and human genetics.
> Forty-four percent are hard-to-drug. Ten are clinical-grade today. IRF1 lights up:
> strong genetics, undruggable. JAK2 and ICOS: green-light. That's a target shortlist
> triaged by *feasibility* — not a ranked gene list.
>
> We surveyed ninety-one prior methods — fifteen years of them. Fourteen use measured
> data. *Zero* return a feasibility verdict. Zero can prove something is out of reach.
> This whole quadrant — measured *and* feasible — was empty. Now it isn't.
>
> And it travels. The same operator, unchanged, jumps from T cells to a completely
> different screen — different cell type, opposite direction — and still recovers a
> held-out target at cosine point-eight-eight. Every number ships with a null model, a
> confidence score, and a one-command reproduce.
>
> Seven interactive explorers, a full manuscript, a ninety-one-method survey — built
> in a week with Claude Science. We answer the question that saves the most money in
> drug discovery: before you run it — can you even get there?

---

## Production checklist

**Highest-leverage move:** screen-record the live HTML explorers, don't use the
static PNGs where you can avoid it. They're dark-themed, animate on hover/scroll, and
read as "a real tool," which is exactly the impression you want.

1. **Capture (in order):** `explorers/pharma_funnel.html` →
   `explorers/reachability_explorer.html` (Panels A, B) →
   `explorers/pharma_triage.html` (hover IRF1/JAK2/ICOS) →
   `explorers/pharma_capability.html` (the empty quadrant) → `index.html`
   (the narrative walkthrough, scroll through all six chapters).
   All open by double-click from `file://` — or serve locally with
   `cd app && python -m http.server 8000` — record the browser at 1920×1080,
   hide bookmarks bar.
2. **Stills for the "paper" beats:** `docs/figures/fig_central_illustration.png`
   (Panels A, B, C, D) — use slow Ken-Burns push-ins, ~4 s each.
3. **Audio:** record VO first, cut picture to match (not the reverse). ElevenLabs or
   a clean phone mic in a closet both work. Music: low, neutral, drop it under VO.
4. **On-screen text:** hard-code the three "punch" numbers as lower-thirds —
   **90% fail** · **0 of 91** · **cosine 0.88** — so they land even on mute.
5. **Title card:** project name + logline + "Built with Claude Science" + the CZI
   dataset credit (Zhu et al. 2025, Marson/Pritchard).

## Number-consistency note (read before you film)

The two Th2→Th1 certificate assets report the null at slightly different values.
These are **reconciled** — the explorer now states both bands on the flagship cell —
so you can show them without contradiction:

| Number | `explorers/reachability_explorer.html` | `fig_central_illustration.png` Panel C |
|---|---|---|
| held-out cosine | 0.446 | 0.448 |
| null z | ≈ 45 | ≈ 24 |
| KKT/Farkas residual | 1e-11 | 1.1×10⁻¹¹ |

Both describe **one result**: the cosine gap (0.446 vs 0.448) is rounding, and the
null-z gap is not an inconsistency but two estimates of the same permutation null —
z scales with the number of shuffles (the atlas uses an 8-shuffle band, z ≈ 45; the
manuscript headlines a denser 60-shuffle null, z ≈ 24), while the cosine is
essentially unchanged. The explorer's flagship cell carries a one-line note saying
exactly this. Filming guidance: the cleanest VO is "**tens of sigma above a
permutation null**" without a specific z; if you show a z, the explorer's inline note
keeps it self-explanatory. The **39% / 25% split and the named activation-certificate
genes** (LYAR, IKZF3, CRTAM…) agree across both assets, and the **cosine 0.88
transfer** (Beat 7) remains your safest headline metric.

## Judge-alignment notes

- **Gladstone Special Prize** ("most potential to advance science that can overcome
  disease"): keep the *autoimmune / T-cell disease* thread explicit — the triage
  shortlist and the disease-genetics cross-reference are your direct hooks. Say
  "autoimmune" out loud at least once.
- **Research-track reproducibility:** Beat 7's "null model + confidence + one-command
  reproduce" line is doing real work — it signals rigor to a scientific panel. Don't cut it.
- **"Built with Claude":** the closing beat credits the workflow (method + 91-method
  survey + manuscript + explorers in a week). On-brand and worth 3 seconds.

## Optional 60-second cutdown

Beats **1 → 4 → 6 → 8** only (hook → verdict+certificate → empty quadrant → close).
That keeps the reframe, the one falsifiable result, the novelty proof, and the CTA.
