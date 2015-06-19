deps = $(wildcard *.h) $(wildcard ../*.h) Makefile

csrcs = $(wildcard *.c)
progs = $(patsubst %.c, %, $(csrcs))

all: $(progs)

$(progs): % : %.c $(deps)
	gcc -Wall -Wextra -O3 $< -lm -o $@

Bossman:
	make -C .. Bossman

lonestar:
	make -C .. lonestar

clean:
	make -C .. clean

