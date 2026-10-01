#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
mkdir -p build
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build livetime.tex > build/pass1.stdout
pdflatex -interaction=nonstopmode -halt-on-error -file-line-error -output-directory=build livetime.tex > build/pass2.stdout
cp build/livetime.pdf KinC_x60_4b_LH2_livetime_Beamer.pdf
pdftotext -layout KinC_x60_4b_LH2_livetime_Beamer.pdf build/slides.txt
printf 'Built KinC_x60_4b_LH2_livetime_Beamer.pdf\n'
