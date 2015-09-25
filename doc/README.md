# Overview

`whamdiis.tex` is the main manuscript, written in LaTeX format.

`whamdiis_notes.tex` is a collection of notes for the manuscript, also in LaTeX format.


## Manuscript related files

File                | Description
--------------------|----------------------------------
whamdiis.tex        | manuscript
whamdiis_notes.tex  | notes for the manuscripts
simul.bib           | reference library
gMOS.bst            | Molecular Simulation style file
review1.txt         | reviewer's report
CHANGES             | changes made in response to the referee's report
/fig                | gray-scale figures
/figclr             | color figures (counterparts of `fig`)
/mathematica        | mathematica scripts to help calculation



## Compiling the manuscript

To automatically compile the manuscript, run
```
make
```
or `make whamdiis.pdf`


Run the following to manually compile the manuscript
```
pdflatex whamdiis
bibtex whamdiis
pdflatex whamdiis
pdflatex whamdiis
```

Change `whamdiis` to `whamdiis_notes` to compile the notes


## Recompile figures

```
make -C fig -B
```


## Make a zip file

```
make zip
```



## Utility scripts

`deannote.py` can be used to remove annotations in the tex file.

`rmbiburl.py` can be used to remove url and doi entries in .bib file
