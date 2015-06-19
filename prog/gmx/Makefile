deps = $(wildcard *.h) $(wildcard ../*.h) Makefile

csrcs = $(wildcard *.c)
progs = $(patsubst %.c, %, $(csrcs))

all: $(progs)

$(progs): % : %.c $(deps)
	gcc -O3 -Wall -Wextra $< -lm -o $@

# use gcc for reproducibility
#icc -Wall -Wremarks $< -o $@

Bossman:
	make -C .. Bossman

lonestar:
	make -C .. lonestar

clean:
	make -C .. clean

