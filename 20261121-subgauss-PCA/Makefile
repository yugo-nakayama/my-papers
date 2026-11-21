PAPER=main
LATEXMK=latexmk

.PHONY: all clean arxiv simulations

all:
	TEXINPUTS=style//: $(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error $(PAPER).tex

clean:
	$(LATEXMK) -C $(PAPER).tex
	rm -f $(PAPER).bbl $(PAPER).run.xml

simulations:
	python3 code/run_simulations.py --out .simulation-output --reps 1000
	cp .simulation-output/simulation_summary.csv results/
	cp .simulation-output/manifest.json results/
	cp .simulation-output/figures/*.pdf figures/

arxiv: all
	bash package_arxiv.sh
