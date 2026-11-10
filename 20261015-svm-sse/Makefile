TARGET = main
ARCHIVE = arxiv_upload

.PHONY: all pdf arxiv sim sim-all clean

all: pdf

pdf: $(TARGET).pdf

$(TARGET).pdf: $(TARGET).tex sections/*.tex figures/*
	pdflatex -interaction=nonstopmode $(TARGET)
	pdflatex -interaction=nonstopmode $(TARGET)

# arXiv提出用にファイルをフラット化してまとめるターゲット
arxiv:
	rm -rf $(ARCHIVE)
	mkdir -p $(ARCHIVE)
	# ソースファイルと分割TeXファイルをコピー
	cp $(TARGET).tex $(ARCHIVE)/
	cp sections/*.tex $(ARCHIVE)/ 2>/dev/null || :
	# references.bib があればコピー
	cp references.bib $(ARCHIVE)/ 2>/dev/null || :
	# cls/styファイルがあればコピー
	cp *.cls *.sty style/*.cls style/*.sty $(ARCHIVE)/ 2>/dev/null || :
	# figuresの中身をフラットにコピー
	cp figures/* $(ARCHIVE)/ 2>/dev/null || :
	# sections/ 参照をフラットな参照（直下）に置換 (main.tex内の \input{sections/...} -> \input{...})
	sed -i '' 's|sections/||g' $(ARCHIVE)/$(TARGET).tex || sed -i 's|sections/||g' $(ARCHIVE)/$(TARGET).tex
	# graphicspath{{figures/}} を graphicspath{{./}} に置換
	sed -i '' 's|{figures/}|{./}|g' $(ARCHIVE)/$(TARGET).tex || sed -i 's|{figures/}|{./}|g' $(ARCHIVE)/$(TARGET).tex
	# tar.gzに圧縮
	cd $(ARCHIVE) && tar -czvf ../arxiv_submission.tar.gz *
	@echo "=== Generated arxiv_submission.tar.gz ==="

sim:
	python3 spiked_svm_sim.py

sim-all: sim
	python3 spiked_svm_gram.py
	python3 spiked_svm_inconsistency.py
	python3 spiked_svm_bias.py
	python3 spiked_svm_classifiers.py

clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.synctex.gz
