#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
mkdir -p build
for doc in report presentation; do
  for pass in 1 2 3; do
    pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build "$doc.tex" > "build/${doc}_pass${pass}.stdout"
  done
done
cp build/report.pdf KinC_x60_4b_LH2_livetime_report.pdf
cp build/presentation.pdf KinC_x60_4b_LH2_livetime_presentation.pdf
pdftotext -layout build/report.pdf build/report.txt
pdftotext -layout build/presentation.pdf build/presentation.txt
