TARGET = main
ARCHIVE = arxiv_upload

.PHONY: all clean arxiv

all:
	latexmk -pdf $(TARGET)

clean:
	latexmk -c

# arXiv提出用にサブディレクトリの図などをフラットにまとめるターゲット
arxiv:
	rm -rf $(ARCHIVE)
	mkdir -p $(ARCHIVE)
	# ソースファイルとbibファイルをコピー
	cp $(TARGET).tex $(ARCHIVE)/
	cp references.bib $(ARCHIVE)/ 2>/dev/null || :
	# cls/sty/bstファイルがあればコピー
	cp *.cls *.sty *.bst style/*.cls style/*.sty style/*.bst $(ARCHIVE)/ 2>/dev/null || :
	# figuresの中身をフラットにコピー（ファイル名衝突に注意）
	cp figures/* $(ARCHIVE)/ 2>/dev/null || :
	# bblファイル（arXivはbbl直納入が安全な場合が多い）を生成してコピー
	latexmk -pdf $(TARGET) || :
	cp $(TARGET).bbl $(ARCHIVE)/ 2>/dev/null || :
	# tar.gzに圧縮
	cd $(ARCHIVE) && tar -czvf ../arxiv_submission.tar.gz *
	@echo "=== Generated arxiv_submission.tar.gz ==="

