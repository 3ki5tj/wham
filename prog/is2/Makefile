deps = $(wildcard *.h) $(wildcard ../*.h) Makefile

csrcs = $(wildcard *.c)
progs = $(patsubst %.c, %, $(csrcs))

all: $(progs)

$(progs): % : %.c $(deps)
	gcc -Wall -Wextra $< -lm -o $@

remote = cz1@129.109.88.204
remoteBossman = $(remote):/Bossman/cz1/mu

Bossman: clean
	rsync -avLz *.[ch] *.py Makefile $(remoteBossman)/prog/ie/

clean:
	find . -name "*~" | xargs rm -f
	rm -rf bin obj a.out *.dat $(progs) \
		*.obj *.exe
	rstrip.py -Rlv

