LATEXMK=latexmk
MAIN=main

.PHONY: all clean sim sim-quick arxiv
all:
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex

sim:
	python3 experiments/run_experiments.py

sim-quick:
	python3 experiments/run_experiments.py --quick

arxiv: all
	rm -rf arxiv_submission arxiv_submission.tar.gz
	mkdir -p arxiv_submission/sections arxiv_submission/figures arxiv_submission/results arxiv_submission/experiments
	cp main.tex main.bbl references.bib arxiv_submission/
	cp sections/*.tex arxiv_submission/sections/
	cp figures/*.pdf arxiv_submission/figures/
	cp results/*_table.tex results/*.csv results/run_metadata.json arxiv_submission/results/
	cp experiments/run_experiments.py experiments/requirements.txt arxiv_submission/experiments/
	cp README.md arxiv_submission/
	tar -czf arxiv_submission.tar.gz -C arxiv_submission .

clean:
	$(LATEXMK) -C
	rm -f main.bbl main.run.xml
