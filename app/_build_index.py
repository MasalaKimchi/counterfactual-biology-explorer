#!/usr/bin/env python3
"""Build app/index.html — the single-page narrative that threads all seven
explorers into one story with interpretive prose. Each explorer is embedded
live via <iframe srcdoc> (fully inlined: works on double-click, no server,
no cross-origin file:// blocking) and is ALSO kept standalone under
app/explorers/ for full-screen viewing. Re-run after editing any explorer:
    python app/_build_index.py
"""
import html, os, pathlib

HERE = pathlib.Path(__file__).resolve().parent
EXP = HERE / "explorers"

def load(name):
    return (EXP / f"{name}.html").read_text(encoding="utf-8")

# ---- chapters: (id, kicker, title, lede, explorer(s), what_you_see, what_it_means) ----
# Each explorer entry is (explorer_id, min_height_px).
CHAPTERS = [
  dict(
    id="problem", num="01", kicker="The problem",
    title="Drug programs fail late, and mostly for the same reason",
    lede=("About nine in ten clinical programmes never reach approval, and the failures cluster "
          "where they cost the most. Phase II is the lowest transition in the whole pipeline — "
          "roughly 30.7% — and about half of late-stage failures are for lack of efficacy: the "
          "target simply was not a way to move the biology. Every AI tool in this space is a "
          "<em>forward</em> predictor — it scores what a perturbation would do. None can tell a "
          "team, from measured data, that a desired cell-state change is out of reach for the "
          "chosen modality. So unreachability is discovered after the screen, or after the trial."),
    explorers=[("pharma_funnel", 620)],
    see=("The clinical attrition funnel, band width proportional to cumulative survival. Hover a "
         "stage to see which of four value levers a feasibility verdict pulls there."),
    mean=("The leak is concentrated in the middle of the funnel, and human genetic support already "
          "carries ~2.6&times; the approval odds (Minikel 2024). A test that flags unreachable "
          "targets <em>before</em> the screen attacks the funnel exactly where it leaks — the "
          "motivation for everything below."),
  ),
  dict(
    id="reframe", num="02", kicker="The reframe",
    title="Ask a different question of data you already have",
    lede=("The fix is not a bigger model. A genome-scale Perturb-seq screen already measures an "
          "average treatment effect for every gene knockdown. Stack those measured effect vectors "
          "and one geometric question falls out: does the target state lie inside the cone of "
          "everything the knockdowns can reach? That is a counterfactual query answered by "
          "design-based causal inference — using constructs the analysis pipeline already builds."),
    explorers=[("causal_reframe", 640)],
    see=("The pipeline as a chain — Perturb-seq &rarr; average treatment effect &rarr; non-negative "
         "cone &rarr; counterfactual query. Click any code object to see the causal-inference name "
         "it already implements."),
    mean=("Nine of ten constructs the reframe needs were already in hand; the reframe &ldquo;costs "
          "nothing.&rdquo; The novelty is the <em>question</em> — feasibility, not prediction — not "
          "a new estimator."),
  ),
  dict(
    id="verdict", num="03", kicker="The instrument",
    title="A falsifiable verdict on real CRISPRi effects",
    lede=("Here is the oracle on the real 33,983-knockdown &times; 10,282-gene matrix. For the "
          "flagship Th2&rarr;Th1 switch in resting CD4&#8314; T cells the verdict is "
          "<em>partially reachable</em>: knockdown can remove the Th2 program but cannot install "
          "the Th1 one. The output is not a similarity score — it is a signed decomposition, a "
          "minimal knockdown recipe, and a Farkas/KKT certificate naming the genes no knockdown "
          "mix can deliver."),
    explorers=[("reachability_explorer", 900)],
    see=("Left: the target's exact shadow in the plane of {reachable fit, residual} — the 2-D "
         "picture is the true projection of the high-dimensional NNLS fit, not a schematic. Right: "
         "pick any of the 12 atlas cells; read its verdict, the signed LOF/GOF/neither split, the "
         "greedy recipe, and the activation certificate."),
    mean=("For the flagship cell: <strong>39% reachable by knockdown, 25% provably needs "
          "activation, 35% neither</strong>; held-out cosine 0.446 clearing its permutation null, "
          "KKT/Farkas residual 1&times;10&#8315;&#185;&#185;. Master-regulator controls land "
          "correctly (GATA3 recovered; TBX21 anti-aligned under knockdown), and across all 12 "
          "atlas cells knockdown is <em>never</em> the majority modality — activation and "
          "irreducible components always dominate. That last fact is the honest headline: the "
          "method is most useful for saying what a knockdown screen <em>cannot</em> do."),
  ),
  dict(
    id="trust", num="04", kicker="The trust layer",
    title="A verdict is only as good as its assumptions — so test them",
    lede=("A counterfactual claim rests on identifying assumptions. Expert causal inference is "
          "making each one explicit, testing the testable ones with data in hand, and quantifying "
          "how far a violation must go before the verdict flips. Two views do this: the six-assumption "
          "stack, and the instrumental-variable treatment of imperfect knockdown."),
    explorers=[("causal_trust", 780), ("causal_iv", 760)],
    see=("First panel: the six identifying assumptions, each with its stress-test and a "
         "robust / calibrated / caveat grade. Second panel: guide assignment as a randomized "
         "instrument and realized knockdown as the treatment — toggle intent-to-treat against the "
         "compliance-rescaled (LATE) verdict."),
    mean=("The headline verdict is robust to the assumptions the design controls: a convex cone is "
          "invariant to per-generator rescaling, so the ITT and LATE verdicts coincide to machine "
          "precision (|&Delta;cosine| = 2&times;10&#8315;&#185;&#8310;), and dropping invalid "
          "instruments moves it by &le;0.0004. The stack is honest about what it does <em>not</em> "
          "yet control — cytokine spillover (SUTVA) is named as the primary external-validity "
          "caveat, and the effect-homogeneity test that needs raw single-cell counts is flagged as "
          "the single most useful next build."),
  ),
  dict(
    id="impact", num="05", kicker="The payoff",
    title="From a verdict to a GO / STOP / REDIRECT shortlist",
    lede=("A feasibility verdict earns its keep when it changes a decision. Take the 102 knockdown "
          "nominations the oracle says are actually needed for the flagship axis and cross them "
          "against Open Targets tractability and immune-disease genetics."),
    explorers=[("pharma_triage", 900)],
    see=("Each point is a required knockdown node, placed by druggability against genetic support. "
         "Filter by triage tier; the three tiers are the GO / STOP / REDIRECT buckets."),
    mean=("45 of 102 (44%) are hard-to-drug and only 10 are clinical-grade today. <strong>IRF1</strong> "
          "is the headline collision — a top-genetics node with no conventional drug handle "
          "(REDIRECT); <strong>JAK2</strong> and <strong>ICOS</strong> are green-lit (GO). This is a "
          "shortlist triaged by feasibility, not a ranked gene list — the deliverable a target team "
          "can act on."),
  ),
  dict(
    id="novelty", num="06", kicker="The delta",
    title="This capability was missing from fifteen years of methods",
    lede=("Is the feasibility verdict actually new? A survey of 91 prior methods (92 including this "
          "work), 2011&ndash;2026, places every one on two axes: does it use measured effects, and "
          "does it return an achievability verdict?"),
    explorers=[("pharma_capability", 760)],
    see=("The capability landscape. 14 methods use measured data; all of them only predict or rank. "
         "The measured &times; achievability quadrant — top-right — is empty."),
    mean=("Zero of 91 return a feasibility verdict on measured effects. This work is the sole "
          "occupant of that quadrant. The novelty is claimed for the method and its decision layer, "
          "not for the regulator biology it recovers — recovering GATA3/TBX21 is validation that "
          "the instrument works, and the same operator transfers unchanged to a K562 CRISPRa screen "
          "(held-out CEBPA at cosine 0.878)."),
  ),
]

FOOT_LINKS = [
  ("The paper (preprint)", "../manuscript/main.pdf"),
  ("Technical Dossier (162 pp)", "../docs/Technical_Dossier.pdf"),
  ("Reproduce", "../reproduce.sh"),
  ("Repository map", "../README.md"),
]

def esc_srcdoc(s):
    # srcdoc attribute value: escape & and " (and <,> harmlessly) — html.escape does all.
    return html.escape(s, quote=True)

def chapter_html(c, idx):
    frames = ""
    for exp_id, minh in c["explorers"]:
        src = esc_srcdoc(load(exp_id))
        title = {
          "reachability_explorer":"Reachability Explorer",
          "causal_reframe":"The reframe that costs nothing",
          "causal_iv":"IV / convex-cone invariance",
          "causal_trust":"The six-assumption trust stack",
          "pharma_funnel":"Attrition funnel + value levers",
          "pharma_triage":"Modality triage — 102 nominations",
          "pharma_capability":"0 of 91 — the empty quadrant",
        }[exp_id]
        frames += f'''
      <div class="frame">
        <div class="framebar">
          <span class="live"><i></i>live &middot; {title}</span>
          <a href="explorers/{exp_id}.html" target="_blank" rel="noopener">open standalone &#8599;</a>
        </div>
        <iframe class="explorer" data-minh="{minh}" loading="lazy"
                srcdoc="{src}" title="{title}"></iframe>
      </div>'''
    return f'''
  <section class="chap" id="{c['id']}">
    <div class="chead">
      <span class="num">{c['num']}</span>
      <div>
        <p class="kicker">{c['kicker']}</p>
        <h2>{c['title']}</h2>
      </div>
    </div>
    <p class="lede">{c['lede']}</p>
    {frames}
    <div class="readout">
      <div class="ro"><h3>What you&rsquo;re looking at</h3><p>{c['see']}</p></div>
      <div class="ro"><h3>What it means</h3><p>{c['mean']}</p></div>
    </div>
  </section>'''

nav = "".join(f'<a href="#{c["id"]}"><b>{c["num"]}</b>{c["kicker"]}</a>' for c in CHAPTERS)
chapters = "".join(chapter_html(c, i) for i, c in enumerate(CHAPTERS))
footlinks = "".join(f'<a href="{href}">{label}</a>' for label, href in FOOT_LINKS)

DOC = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>Cell-State Reachability — an interactive walkthrough</title>
<style>
  :root{{
    --bg:#0f1419; --panel:#161c24; --ink:#e8edf2; --dim:#8b97a6; --line:#2a333f;
    --lof:#1f9e89; --gof:#d99b2b; --neither:#c7cdd6; --accent:#4aa3df; --good:#4fd6c0;
  }}
  *{{box-sizing:border-box}}
  html{{scroll-behavior:smooth}}
  body{{margin:0;background:var(--bg);color:var(--ink);
    font:15px/1.62 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}}
  a{{color:var(--accent);text-decoration:none}} a:hover{{text-decoration:underline}}
  em{{font-style:italic;color:#c7d2de}} strong{{color:#fff;font-weight:650}}
  code{{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;color:#9fb0c0;font-size:.92em}}

  /* hero */
  .hero{{max-width:1100px;margin:0 auto;padding:64px 26px 30px}}
  .tag{{font-size:11px;text-transform:uppercase;letter-spacing:1.6px;color:var(--accent);margin:0 0 14px;font-weight:650}}
  .hero h1{{font-size:34px;line-height:1.16;font-weight:680;letter-spacing:.2px;margin:0 0 16px;max-width:20ch}}
  .hero p.big{{font-size:17px;line-height:1.6;color:#c7d2de;max-width:70ch;margin:0 0 14px}}
  .hero p.small{{font-size:13.5px;color:var(--dim);max-width:72ch;margin:0}}
  .hero .kpis{{display:flex;gap:26px;flex-wrap:wrap;margin:26px 0 0}}
  .hero .kpi b{{display:block;font-size:26px;font-weight:680;font-variant-numeric:tabular-nums;letter-spacing:.3px}}
  .hero .kpi span{{font-size:11.5px;color:var(--dim);text-transform:uppercase;letter-spacing:.7px}}

  /* sticky nav */
  nav{{position:sticky;top:0;z-index:20;background:rgba(15,20,25,.86);backdrop-filter:blur(9px);
    border-top:1px solid var(--line);border-bottom:1px solid var(--line)}}
  nav .in{{max-width:1100px;margin:0 auto;padding:9px 22px;display:flex;gap:4px;flex-wrap:wrap}}
  nav a{{display:flex;align-items:baseline;gap:7px;color:var(--dim);font-size:12.5px;
    padding:6px 11px;border-radius:7px;white-space:nowrap}}
  nav a b{{color:var(--accent);font-variant-numeric:tabular-nums;font-size:11px}}
  nav a:hover{{color:var(--ink);background:var(--panel);text-decoration:none}}
  nav a.on{{color:var(--ink);background:#182533}}

  /* chapters */
  .chap{{max-width:1100px;margin:0 auto;padding:52px 26px 10px;border-top:1px solid var(--line)}}
  .chap:first-of-type{{border-top:none}}
  .chead{{display:flex;gap:18px;align-items:flex-start;margin:0 0 14px}}
  .num{{font-size:34px;font-weight:700;color:var(--line);font-variant-numeric:tabular-nums;line-height:1;flex:none}}
  .kicker{{font-size:11px;text-transform:uppercase;letter-spacing:1.5px;color:var(--accent);margin:2px 0 3px;font-weight:650}}
  .chap h2{{font-size:23px;font-weight:660;letter-spacing:.2px;margin:0;line-height:1.25}}
  .lede{{font-size:15.5px;line-height:1.7;color:#c2cdd8;max-width:78ch;margin:0 0 24px}}

  /* embedded explorer */
  .frame{{border:1px solid var(--line);border-radius:11px;overflow:hidden;margin:0 0 18px;background:#0b0f14;
    box-shadow:0 10px 34px rgba(0,0,0,.32)}}
  .framebar{{display:flex;justify-content:space-between;align-items:center;padding:8px 14px;
    background:linear-gradient(180deg,#12181f,#0e131a);border-bottom:1px solid var(--line);font-size:12px}}
  .framebar .live{{color:var(--dim);letter-spacing:.4px;display:flex;align-items:center;gap:7px}}
  .framebar .live i{{width:7px;height:7px;border-radius:50%;background:var(--good);display:inline-block;
    box-shadow:0 0 0 0 rgba(79,214,192,.6);animation:pulse 2.4s infinite}}
  @keyframes pulse{{0%{{box-shadow:0 0 0 0 rgba(79,214,192,.5)}}70%{{box-shadow:0 0 0 7px rgba(79,214,192,0)}}100%{{box-shadow:0 0 0 0 rgba(79,214,192,0)}}}}
  iframe.explorer{{width:100%;border:0;display:block;background:var(--bg)}}

  /* readout */
  .readout{{display:grid;grid-template-columns:1fr 1fr;gap:22px;margin:0 0 8px;padding:2px 0 6px}}
  @media(max-width:760px){{.readout{{grid-template-columns:1fr}}}}
  .ro h3{{font-size:11px;text-transform:uppercase;letter-spacing:1.2px;color:var(--dim);margin:0 0 6px;font-weight:600}}
  .ro p{{font-size:14px;line-height:1.66;color:#b9c4d0;margin:0}}
  .ro:first-child{{border-left:2px solid var(--accent);padding-left:16px}}
  .ro:last-child{{border-left:2px solid var(--lof);padding-left:16px}}

  footer{{max-width:1100px;margin:34px auto 0;padding:26px;border-top:1px solid var(--line);color:var(--dim);font-size:12.5px}}
  footer .fl{{display:flex;gap:20px;flex-wrap:wrap;margin:0 0 14px}}
  footer .fl a{{font-size:13.5px}}
  footer p{{margin:0;max-width:80ch;line-height:1.6}}
</style>
</head>
<body>

<div class="hero">
  <p class="tag">Cell-State Reachability &middot; Built with Claude &mdash; Life Sciences</p>
  <h1>What a knockdown screen <em>can&rsquo;t</em> reach &mdash; decided before you run it.</h1>
  <p class="big">Every AI model in biology predicts what a perturbation <em>will</em> do. This one
  returns a falsifiable <strong>reachable / provably-out-of-reach</strong> verdict &mdash; with the
  minimal knockdown recipe and a numerical certificate &mdash; for a real genome-scale CRISPRi
  Perturb-seq screen in primary human CD4&#8314; T cells.</p>
  <p class="small">Below is the whole argument as one walkthrough: why it matters, the reframe that
  makes it work, the instrument on real data, why you can trust the verdict, what it changes for a
  drug program, and why the capability was missing. Each panel is the live explorer &mdash; interact
  with it in place, or open it full-screen.</p>
  <div class="kpis">
    <div class="kpi"><b style="color:var(--lof)">39/25/35</b><span>LOF / GOF / neither</span></div>
    <div class="kpi"><b>0.446</b><span>held-out cosine</span></div>
    <div class="kpi"><b style="color:var(--gof)">0 of 91</b><span>priors with a verdict</span></div>
    <div class="kpi"><b>102</b><span>nominations triaged</span></div>
  </div>
</div>

<nav><div class="in">{nav}</div></nav>

{chapters}

<footer>
  <div class="fl">{footlinks}</div>
  <p>Every number in every panel is computed by <code>reachability.py</code> on the real
  <code>GWCD4i.DE_stats.h5ad</code> effect matrix and cross-checked against the committed
  <code>results/</code> tables; each explorer embeds its own verified data payload and states its
  provenance in its own footer. Dataset: genome-scale CRISPRi Perturb-seq in primary human CD4&#8314;
  T cells (Zhu et al. 2025, Marson &amp; Pritchard labs, CZI Virtual Cells Platform). The seven
  explorers also live standalone under <code>app/explorers/</code>.</p>
</footer>

<script>
// Auto-size each srcdoc iframe to its content (same-origin, so contentDocument is readable),
// and keep it sized as the embedded explorer reflows or the user interacts.
function fit(f){{
  try{{
    var d = f.contentDocument || f.contentWindow.document;
    var h = Math.max(d.documentElement.scrollHeight, d.body ? d.body.scrollHeight : 0);
    var min = parseInt(f.getAttribute('data-minh')||'0',10);
    f.style.height = Math.max(h, min) + 'px';
  }}catch(e){{ f.style.height = (f.getAttribute('data-minh')||600) + 'px'; }}
}}
document.querySelectorAll('iframe.explorer').forEach(function(f){{
  f.addEventListener('load', function(){{
    fit(f);
    try{{
      var d = f.contentDocument || f.contentWindow.document;
      if (window.ResizeObserver && d.body){{
        var ro = new ResizeObserver(function(){{ fit(f); }});
        ro.observe(d.body);
      }}
      // atlas/tab clicks inside the explorer change height — remeasure shortly after
      d.addEventListener('click', function(){{ setTimeout(function(){{ fit(f); }}, 380); }}, true);
    }}catch(e){{}}
  }});
  // initial guess before load
  f.style.height = (f.getAttribute('data-minh')||600) + 'px';
}});
window.addEventListener('resize', function(){{
  document.querySelectorAll('iframe.explorer').forEach(fit);
}});

// nav active-state on scroll
var secs = [].slice.call(document.querySelectorAll('.chap'));
var links = {{}};
document.querySelectorAll('nav a').forEach(function(a){{ links[a.getAttribute('href').slice(1)] = a; }});
if ('IntersectionObserver' in window){{
  var io = new IntersectionObserver(function(es){{
    es.forEach(function(e){{
      if(e.isIntersecting){{
        Object.values(links).forEach(function(a){{a.classList.remove('on');}});
        if(links[e.target.id]) links[e.target.id].classList.add('on');
      }}
    }});
  }}, {{rootMargin:'-45% 0px -50% 0px'}});
  secs.forEach(function(s){{ io.observe(s); }});
}}
</script>
</body>
</html>
'''

out = HERE / "index.html"
out.write_text(DOC, encoding="utf-8")
print("wrote", out, os.path.getsize(out), "bytes")
print("chapters:", len(CHAPTERS), "| embedded iframes:", sum(len(c["explorers"]) for c in CHAPTERS))
