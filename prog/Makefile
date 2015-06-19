remote = cz1@129.109.88.204
remoteBossman = $(remote):/Bossman/cz1/wham

Bossman: clean
	rsync -avLz *.[ch]      Makefile        $(remoteBossman)/prog/
	rsync -avLz is2/*.[ch]  is2/Makefile    $(remoteBossman)/prog/is2/
	rsync -avLz lj/*.[ch]   lj/Makefile     $(remoteBossman)/prog/lj/
	rsync -avLz gmx/*.[ch]  gmx/Makefile    $(remoteBossman)/prog/gmx/

remotelonestar = oo1@lonestar.tacc.utexas.edu:scratch/wham

lonestar: clean
	rsync -avLz *.[ch]      Makefile        $(remotelonestar)/prog/
	rsync -avLz is2/*.[ch]  is2/Makefile    $(remotelonestar)/prog/is2/
	rsync -avLz lj/*.[ch]   lj/Makefile     $(remotelonestar)/prog/lj/
	rsync -avLz gmx/*.[ch]  gmx/Makefile    $(remotelonestar)/prog/gmx/

lstar: lonestar

clean:
	find . -name "*~" | xargs rm -f
	find . -name ".*.un~" | xargs rm -f
	rm -rf */bin */obj */a.out */*.dat */*.obj */*.exe
	rstrip.py -Rlv

