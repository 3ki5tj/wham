remote = cz1@129.109.88.204
remoteBossman = $(remote):/Bossman/cz1/wham

Bossman: clean
	rsync -avLz *.py                *.pbs       Makefile        $(remoteBossman)/data/
	rsync -avLz is2/*.py  is2/*.sh              is2/Makefile    $(remoteBossman)/data/is2/
	rsync -avLz lj/*.py                         lj/Makefile     $(remoteBossman)/data/lj/
	rsync -avLz 1VII/*.py 1VII/*.sh 1VII/*.pbs  1VII/Makefile   $(remoteBossman)/data/1VII/

remotelonestar = oo1@lonestar.tacc.utexas.edu:scratch/wham

lonestar: clean
	rsync -avLz *.py                *.pbs       Makefile        $(remotelonestar)/data/
	rsync -avLz is2/*.py  is2/*.sh              is2/Makefile    $(remotelonestar)/data/is2/
	rsync -avLz lj/*.py                         lj/Makefile     $(remotelonestar)/data/lj/
	rsync -avLz 1VII/*.py 1VII/*.sh 1VII/*.pbs  1VII/Makefile   $(remotelonestar)/data/1VII/
	rsync -avLz 1VII/e.ls 1VII/dotdote.ls                       $(remotelonestar)/data/1VII/
	make -C ../prog lonestar

lstar: lonestar

clean:
	find . -name "*~" | xargs rm -f
	find . -name ".*.un~" | xargs rm -f
	rstrip.py -Rlv

