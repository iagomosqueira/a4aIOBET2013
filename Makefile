PDFS:=$(patsubst %.Rmd,%.pdf,$(wildcard *.Rmd))
R_HOME=/usr/local

all: $(PDFS)

clean:
	rm -rf *.tex *.bbl *.blg *.aux *.out *.log *.spl *.md figure/ cache/

cleanall: clean
	rm $(PDFS)

%.pdf: %.md
	pandoc -s -S --biblio $*.bib -o $*.pdf $*.md

%.md: %.Rmd
	Rscript -e "library(knitr)" -e "knit('$*.Rmd')"
