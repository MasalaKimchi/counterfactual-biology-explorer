# Deploying the interactive walkthrough

`index.html` in this folder is the **single-page narrative** — a guided, six-chapter walkthrough
(Problem → Reframe → Verdict → Trust → Impact → Novelty) that embeds all seven interactive explorers
live and threads them with interpretive prose. Each explorer is inlined with an `<iframe srcdoc>`, so
the page is one fully self-contained HTML file: no build step, no server, no external dependencies.

The seven explorers also live standalone in `explorers/` and open independently — the narrative links
to each ("open standalone ↗"), and they are what you screen-record for the demo video.

```
app/
├── index.html          ← the narrative hub (open this)
├── explorers/          ← the 7 standalone explorers
│   ├── reachability_explorer.html
│   ├── causal_reframe.html
│   ├── causal_iv.html
│   ├── causal_trust.html
│   ├── pharma_funnel.html
│   ├── pharma_triage.html
│   └── pharma_capability.html
├── previews/           ← static PNG thumbnails
└── _build_index.py     ← regenerates index.html from the explorers
```

## Fastest path: open it locally

You do not need to deploy anything to use it. Just open the file in any modern browser:

```bash
open app/index.html          # macOS
xdg-open app/index.html      # Linux
start app\index.html         # Windows
```

It works offline, from a `file://` URL — no internet connection required. Because each explorer is
embedded as an inline `srcdoc` (not a linked `src`), the iframes render live even from `file://`,
with no cross-origin restriction and no local web server.

## Regenerating the narrative

`index.html` is generated from the seven explorers and the interpretive prose in `_build_index.py`.
After editing any explorer (or the prose), rebuild it:

```bash
python app/_build_index.py     # reads explorers/*.html, rewrites app/index.html
```

## Publish on GitHub Pages (free static hosting)

Because everything is static files, GitHub Pages serves them as-is. Naming the hub `index.html` means
the folder's URL resolves straight to the walkthrough.

1. **Push the repo** (this `app/` folder included) to GitHub, on the branch you want to publish from
   (usually `main`).
2. **Enable Pages.** In the repository on github.com, go to **Settings → Pages**.
   - Under **Build and deployment → Source**, choose **Deploy from a branch**.
   - Set **Branch** to `main` and the **folder** to **`/app`**, then **Save.**
     *(If your account only offers `/ (root)` and `/docs`, either copy the `app/` contents to a
     top-level `docs/` folder and pick `/docs`, or publish from the repo root and pick `/ (root)`.)*
3. **Wait ~1–2 minutes** for the first build. The Pages panel then shows a green check and the live URL.

### Resulting URL

```
https://<your-username>.github.io/<your-repo-name>/
```

For example, if the GitHub user is `octocat` and the repo is `cell-state-reachability`:

```
https://octocat.github.io/cell-state-reachability/
```

Because the file is named `index.html`, that bare folder URL loads the explorer directly. Every push to
the publish branch re-deploys automatically.

## Why no build step

| Concern | Status |
|---|---|
| External JS/CSS/CDN dependencies | **None** — no `<script src>`, no `<link href>`, no `@import`. |
| Data fetch at runtime | **None** — the dataset is embedded inline as `const DATA`; no `fetch`/XHR. |
| Web fonts / external images | **None** — system fonts and inline SVG only. |
| Server / API needed | **No** — it is a static document. |

You can confirm self-containment yourself:

```bash
# should print nothing (no external references)
grep -Ei '<script[^>]*src=|<link[^>]*href=|https?://|@import|fetch\(' deploy/index.html
```

*(The only `url(...)` occurrences in the file are internal SVG fragment references like `url(#ah)` —
these point at `<defs>` inside the same document and are offline-safe.)*

## Custom domain (optional)

To serve it from your own domain, add a `CNAME` file next to `index.html` containing the domain (e.g.
`explorer.example.org`) and configure the DNS record per GitHub's
[custom-domain docs](https://docs.github.com/pages/configuring-a-custom-domain-for-your-github-pages-site).
