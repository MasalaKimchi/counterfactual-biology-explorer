#!/usr/bin/env bash
# Compile the CombiCone manuscript to PDF.
# Requires a working TeX installation (TeX Live / MacTeX / Overleaf).
# One command:  bash compile.sh
set -euo pipefail
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode main.tex
bibtex main || true
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
echo "Built main.pdf"
