deps = $(wildcard *.h) $(wildcard ../*.h) Makefile

csrcs = $(wildcard *.c)
progs = $(patsubst %.c, %, $(csrcs))
progs_ldbl = $(patsubst %, %_ldbl, $(filter-out xvgact, $(progs)))


all: $(progs) $(progs_ldbl)

$(progs): % : %.c $(deps)
	gcc -O3 -Wall -Wextra $< -lm -o $@

$(progs_ldbl): %_ldbl : %.c $(deps)
	gcc -DLDBL -O3 -Wall -Wextra $< -lm -o $@

# use gcc for reproducibility
#icc -Wall -Wremarks $< -o $@

Bossman:
	make -C .. Bossman

lonestar:
	make -C .. lonestar

clean:
	$(RM) $(progs) $(progs_ldbl)
	make -C .. clean

