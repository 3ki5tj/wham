prj=whamdiis

figgps=$(wildcard fig/*.gp)
pdffigs=$(patsubst %.gp,%.pdf, $(figgps))
epsfigs=$(patsubst %.gp,%.eps, $(figgps))

texpdf = $(prj).pdf $(prj)_notes.pdf

all: $(texpdf)

$(texpdf) : %.pdf : %.tex simul.bib $(pdffigs) $(figgps) Makefile
	make -C fig
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi \
	      *Notes.bib *~
	rstrip.py

$(prj)doc.zip: $(prj).tex $(prj).pdf \
  	       $(prj)_notes.tex $(prj)_notes.pdf \
  	       simul.bib \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       Makefile
	zip -r $@ $^

zip: $(prj)doc.zip

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/$(prj)
