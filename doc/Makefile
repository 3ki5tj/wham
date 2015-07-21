prj=whamdiis

ver=3

figgps=$(wildcard fig/*.gp)
pdffigs=$(patsubst %.gp,%.pdf, $(figgps))
epsfigs=$(patsubst %.gp,%.eps, $(figgps))

texpdf = $(prj).pdf $(prj)_notes.pdf

alltexpdf = $(patsubst %.tex,%.pdf, $(wildcard $(prj)*.tex))

regular: $(texpdf)

all: $(alltexpdf)

$(alltexpdf) : %.pdf : %.tex simul.bib $(pdffigs) $(figgps) Makefile
	make -C fig
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi \
	      *Notes.bib *~
	rstrip.py

zip: $(prj)docV$(ver).zip

# make a zip for the current version

$(prj)docV$(ver).zip:
	mkdir -p old/V$(ver)
	rsync -avz $(prj).tex old/V$(ver)/$(prj)_annotatedV$(ver).tex
	rsync -avz $(prj)_notes.tex old/V$(ver)/$(prj)_notesV$(ver).tex
	rsync -avz simul.bib old/V$(ver)/
	rsync -avz Makefile old/V$(ver)/
	rsync -avz README.md old/V$(ver)/
	rsync -avz deannote.py old/V$(ver)/
	cd old/V$(ver) && \
	  ln -sf ../../fig
	cd old/V$(ver) && \
	  python deannote.py -v -i $(prj)_annotatedV$(ver).tex -o $(prj)V$(ver).tex
	make all -C old/V$(ver)/
	cd old/V$(ver) && \
	  zip -r ../../$(prj)V$(ver).zip *.tex *.pdf simul.bib \
	       $(epsfigs) $(pdffigs) $(figgps) Makefile README.md


$(prj)doc.zip: $(prj).tex $(prj).pdf \
  	       $(prj)_notes.tex $(prj)_notes.pdf \
  	       simul.bib \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       Makefile README.md
	zip -r $@ $^

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/$(prj)
