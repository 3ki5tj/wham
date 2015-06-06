deps = $(wildcard *.h) $(wildcard ../*.h) Makefile

csrcs = $(wildcard *.c)
progs = $(patsubst %.c, %, $(csrcs))

all: $(progs)

$(progs): % : %.c $(deps)
	gcc -Wall -O2 $< -lm -o $@

# use gcc for reproducibility
#icc -Wall -Wremarks $< -o $@

remote = cz1@129.109.88.204
remoteBossman = $(remote):/Bossman/cz1/wham

Bossman: clean
	rsync -avLz *.[ch] Makefile $(remoteBossman)/prog/gmx/
	rsync -avLz ../*.[ch] $(remoteBossman)/prog/

remotelonestar = oo1@lonestar.tacc.utexas.edu:scratch/wham

lonestar: clean
	rsync -avLz *.[ch] Makefile $(remotelonestar)/prog/gmx/
	rsync -avLz ../*.[ch] $(remotelonestar)/prog/

clean:
	find . -name "*~" | xargs rm -f
	rm -rf bin obj a.out *.dat $(progs) \
		*.obj *.exe
	rstrip.py -Rlv

