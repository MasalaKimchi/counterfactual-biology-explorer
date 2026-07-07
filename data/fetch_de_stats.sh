#!/usr/bin/env bash
# Fetch the gene-level effect matrix GWCD4i.DE_stats.h5ad (Tier-2).
#
# VERIFIED 2026-07-08: the artifact IS available and directly downloadable from a
# PUBLIC S3 bucket — no CZI login, no vcp-cli, no --file filter needed. The old note
# in this script ("DE_stats is NOT a search hit, cannot cherry-pick by name") is
# obsolete: the dataset is now published as a single collection whose derived
# artifacts sit next to the raw per-donor files in one public prefix.
#
#   Dataset (VCP):   Primary Human CD4+ T Cell Perturb-seq  (v1.0)
#   Dataset ID:      019b479b-e1e8-7c6f-ba55-b6b6936c7a7b
#   Public prefix:   s3://genome-scale-tcell-perturb-seq/marson2025_data/
#   download_access: public   (list + get work with an UNSIGNED / --no-sign-request client)
#
#   File:            GWCD4i.DE_stats.h5ad
#   Size:            16.79 GB  (16,786,240,107 bytes)  <- the FULL h5ad with all layers,
#                    not 1.4 GB. 1.4 GB is one float32 layer once loaded in memory.
#   Grain:           33,983 perturbation x condition rows x 10,282 genes
#   Layers:          log_fc, zscore, p_value, adj_p_value, baseMean, lfcSE
#
# You do NOT need this to run the graded pipeline — Tier-1 (the no-auth CSVs already in
# data/) produces directional nominations on its own. This file is only for the full
# gene-space reachability cone solver (src/reachability.py).
#
# ---------------------------------------------------------------------------
# Method A (RECOMMENDED — no auth, no vcp-cli): direct public-S3 download via boto3.
#   boto3 ships with vcp-cli[data]; if missing: pip install boto3
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$(dirname "$0")/.."

BUCKET="genome-scale-tcell-perturb-seq"
KEY="marson2025_data/GWCD4i.DE_stats.h5ad"
OUT="data/GWCD4i.DE_stats.h5ad"

echo "Downloading s3://${BUCKET}/${KEY}  (~16.8 GB, public, no login) -> ${OUT}"
python3 - "$BUCKET" "$KEY" "$OUT" <<'PY'
import sys, boto3
from botocore import UNSIGNED
from botocore.client import Config
from boto3.s3.transfer import TransferConfig
bucket, key, out = sys.argv[1], sys.argv[2], sys.argv[3]
s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))
size = s3.head_object(Bucket=bucket, Key=key)["ContentLength"]
seen = {"n": 0}
def cb(n):
    seen["n"] += n
    pct = 100 * seen["n"] / size
    print(f"\r  {seen['n']/1e9:6.2f} / {size/1e9:.2f} GB  ({pct:5.1f}%)", end="", flush=True)
s3.download_file(bucket, key, out,
                 Config=TransferConfig(multipart_threshold=64*1024*1024,
                                       max_concurrency=8),
                 Callback=cb)
print("\nDone.")
PY

# ---------------------------------------------------------------------------
# Method B (alternative): CZI vcp-cli. NOTE — as of 0.54.x BOTH `search` and
# `download` now require `vcp login` (a free browser OAuth flow; it did not used to
# gate `search`). vcp still has no per-file-NAME filter, but it DOES now have
# `--extension`, so you can restrict a dataset download to h5ad files:
#
#     pip install 'vcp-cli[data]'
#     vcp login
#     # WARNING: --extension h5ad ALSO matches the 12 raw per-donor *.assigned_guide.h5ad
#     # files (~140-170 GB EACH, ~1.8 TB total). Prefer Method A for just DE_stats.
#     vcp data download --id 019b479b-e1e8-7c6f-ba55-b6b6936c7a7b --extension h5ad -o ./data
#
# ---------------------------------------------------------------------------
# Verify after either method:
#     python -m src.data_loader --check      # -> "DE_stats (h5ad, Tier 2) present: True"
echo
echo "Verify with: python -m src.data_loader --check"
