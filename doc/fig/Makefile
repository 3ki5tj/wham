GP = $(wildcard *.gp)
PNGGP =
EPSGP = $(filter-out $(PNGGP), $(GP))
EPS = $(patsubst %.gp,%.eps, $(EPSGP))
PDF = $(patsubst %.gp,%.pdf, $(EPSGP))
PNG = $(patsubst %.gp,%.png, $(PNGGP))
epsdeps = Makefile $(wildcard ../data/*/*.dat)
pngdeps = Makefile

all: $(PNG) $(PDF)

$(PNG) : %.png : %.gp $(pngdeps)
	gnuplot $<

$(PDF) : %.pdf : %.eps $(epsdeps)
	epstopdf $<

$(EPS) : %.eps : %.gp $(epsdeps)
	gnuplot $<

clean:
	rm -f $(PDF) $(EPS) $(PNG) *~
	rstrip.py -Rv

