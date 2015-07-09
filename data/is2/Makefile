all: mkrundat mktracedat

mkrundat:
	cd run    && ../stat.py && cd ..
	cd runhp  && ../stat.py && cd ..
	cd runhpl && ../stat.py && cd ..
	cd runkth && ../stat.py && cd ..

mktracedat:
	cd trace    && ../stat.py && cd ..
	cd tracehp  && ../stat.py && cd ..
	cd tracehpl && ../stat.py && cd ..
	cd tracekth && ../stat.py && cd ..

lonestar:
	make -C .. lonestar
