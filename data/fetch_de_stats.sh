#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FILE="${1:-GWCD4i.DE_stats.h5ad}"
S3_PREFIX="https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/marson2025_data"
AUTHOR_PREFIX="https://raw.githubusercontent.com/emdann/GWT_perturbseq_analysis_2025/848d62fc2b7027f7218d6fc5f5b0c37255dc94af/metadata/suppl_tables"
OUT="$ROOT/data/$FILE"
PART="$OUT.part"

EXPECTED_SHA256=""
case "$FILE" in
  GWCD4i.DE_stats.h5ad) EXPECTED_BYTES=16786240107; EXPECTED_SHA256=c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62; URL="$S3_PREFIX/$FILE" ;;
  GWCD4i.DE_stats.by_donors.h5mu) EXPECTED_BYTES=16866278447; EXPECTED_SHA256=2ee3cf90925600eb044619021da2bdd47d661f306a204586652256facf17af64; URL="$S3_PREFIX/$FILE" ;;
  GWCD4i.DE_stats.by_guide.h5mu) EXPECTED_BYTES=29424424894; URL="$S3_PREFIX/$FILE" ;;
  GWCD4i.pseudobulk_merged.h5ad) EXPECTED_BYTES=44566657140; URL="$S3_PREFIX/$FILE" ;;
  Th2_Th1_polarization_signature_DE_results_full.suppl_table.csv) EXPECTED_BYTES=6155771; EXPECTED_SHA256=c47d2df21414ca85e7aa255f4148904eec700fbcd9debc2f734ec97049698444; URL="$AUTHOR_PREFIX/$FILE" ;;
  IL10IL21bulkRNAseq_DESeq2_results.csv) EXPECTED_BYTES=13952871; EXPECTED_SHA256=c20418a9285b10104dbae362b825971f86f97425800a92269e4433ce780e666d; URL="$AUTHOR_PREFIX/$FILE" ;;
  IL10_IL21_arrayed_validation.csv) EXPECTED_BYTES=2200; EXPECTED_SHA256=f60cdda392d6f29d10a539727ff7324b04d17e35c0512c889b733e00380b83dc; URL="$AUTHOR_PREFIX/$FILE" ;;
  *)
    echo "Unsupported or unresolved object: $FILE" >&2
    exit 2
    ;;
esac

if [[ -e "$OUT" ]]; then
  echo "Refusing to overwrite existing file: $OUT" >&2
  shasum -a 256 "$OUT"
  exit 2
fi
command -v curl >/dev/null || { echo "curl is required" >&2; exit 2; }
curl --fail --location --continue-at - --output "$PART" "$URL"

ACTUAL_BYTES="$(wc -c < "$PART" | tr -d ' ')"
if [[ "$ACTUAL_BYTES" != "$EXPECTED_BYTES" ]]; then
  echo "Byte-length mismatch: expected $EXPECTED_BYTES, found $ACTUAL_BYTES" >&2
  exit 1
fi
ACTUAL_SHA256="$(shasum -a 256 "$PART" | awk '{print $1}')"
if [[ -n "$EXPECTED_SHA256" && "$ACTUAL_SHA256" != "$EXPECTED_SHA256" ]]; then
  echo "SHA-256 mismatch: expected $EXPECTED_SHA256, found $ACTUAL_SHA256" >&2
  exit 1
fi
mv "$PART" "$OUT"
echo "$ACTUAL_SHA256  $OUT"
