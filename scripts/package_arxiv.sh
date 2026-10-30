#!/usr/bin/env bash
set -euo pipefail

# arXiv提出用にファイルをフラット化してtar.gzにまとめるスクリプト
TARGET="${1:-main}"
ARCHIVE_DIR="arxiv_upload"
OUTPUT_TAR="arxiv_submission.tar.gz"

echo "=== Packaging arXiv submission for target: ${TARGET} ==="

rm -rf "${ARCHIVE_DIR}" "${OUTPUT_TAR}"
mkdir -p "${ARCHIVE_DIR}"

# 1. コンパイルして bbl を生成
latexmk -pdf "${TARGET}" || true

# 2. メインの tex と bib をコピー
cp "${TARGET}.tex" "${ARCHIVE_DIR}/"
if [ -f "references.bib" ]; then
  cp references.bib "${ARCHIVE_DIR}/"
fi

# 3. bbl ファイルをコピー (arXiv は bbl 直納が安全)
if [ -f "${TARGET}.bbl" ]; then
  cp "${TARGET}.bbl" "${ARCHIVE_DIR}/"
fi

# 4. スタイルファイル・クラスファイル (ルート直下および style/ ディレクトリ)
shopt -s nullglob
for f in *.cls *.sty *.bst style/*.cls style/*.sty style/*.bst; do
  [ -f "$f" ] && cp "$f" "${ARCHIVE_DIR}/"
done

# 5. 図版ファイルをフラットにコピー
for fig in figures/*; do
  [ -f "$fig" ] && cp "$fig" "${ARCHIVE_DIR}/"
done
shopt -u nullglob

# 6. tar.gz アーカイブを作成
tar -C "${ARCHIVE_DIR}" -czvf "${OUTPUT_TAR}" .

echo "=== Generated ${OUTPUT_TAR} successfully! ==="
