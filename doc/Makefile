prj=wham

figs=
pdffigs=$(patsubst %,fig/%.pdf, $(figs))
epsfigs=$(patsubst %,fig/%.eps, $(figs))
figgps=$(patsubst %,fig/%.gp,  $(figs))

texpdf = $(prj).pdf

$(texpdf) : %.pdf : %.tex simul.bib Makefile
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi \
	      *Notes.bib *~
	rstrip.py

$(prj)doc.zip: $(prj).tex $(prj).pdf \
  	       simul.bib \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       README Makefile
	zip -r $@ $^

zip: $(prj)doc.zip

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/$(prj)
