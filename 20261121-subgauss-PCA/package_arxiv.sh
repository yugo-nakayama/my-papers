#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
stage="$root/arxiv_submission"
archive="$root/arxiv_submission.tar.gz"
rm -rf "$stage" "$archive"
mkdir -p "$stage/sections" "$stage/figures"
cp "$root/main.tex" "$root/references.bib" "$stage/"
cp "$root"/sections/*.tex "$stage/sections/"
cp "$root"/figures/*.pdf "$stage/figures/"
# Include the generated bibliography so arXiv can compile even without BibTeX reruns.
if [[ -f "$root/main.bbl" ]]; then cp "$root/main.bbl" "$stage/"; fi
(
  cd "$stage"
  TEXINPUTS="$root/style//:" pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
  TEXINPUTS="$root/style//:" pdflatex -interaction=nonstopmode -halt-on-error main.tex >/dev/null
  if grep -Eq 'undefined references|Citation.*undefined' main.log; then
    echo "arXiv staging compile has unresolved references" >&2
    exit 1
  fi
  rm -f main.aux main.log main.out main.pdf main.toc
)
tar -C "$stage" -czf "$archive" .
echo "Created $archive"
