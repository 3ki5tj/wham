Making the comparison figure
=============================

./mkwhamcmp.sh


Protocols
==========

The file give detailed protocols for data collection.

`prep.py` is the script to prepare system for GROMACS.


## NVT ensemble ##

### Prepare systems ###

On lonestar first go to the scratch directory:
```
cds wham/data/1VII
```

Prepare systems one-by-one:
```
./prep.py -T 300 --nt=4
./prep.py -T 310 --nt=4
...
./prep.py -T 410 --nt=4
```
Here, `--nt=4` limits the number of threads to 4.

Locally, the above steps can be combined as
```
./mkT.sh
```



### Run simulations ###

On lonestar
```
qsub T_lonestar.pbs
```

Locally
```
cd T300
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx mdrun -v -deffnm nvt
```



### Export energy ###

```
./T_scan.sh
```

Manually
```
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx energy -f nvt.edr -o e.xvg
```
Select "Potential" (11) and possibly "Pressure" (17).



### Run WHAM ###

```
cd ../../prog/gmx
ls --color=none -d ../../data/1VII/T[0-9][0-9][0-9]/e.xvg > e.ls
./xvgwham e.ls --wham=mdiis
```



### Run MBAR ###

```
./xvgmbar e.ls --wham=mdiis
```



### Benchmark ###

For the number of iteractions needed to reach an error tolerance
```
python xvgrun.py
python stat.py
```

For the error versus the number of iteractions
```
python xvgtrace.py
python stat.py
```





## NPT ensemble ##

### Prepare systems ###

On lonestar first go to the scratch directory:
```
cds wham/data/1VII
```

Prepare systems one-by-one:
```
./prep.py -T 300 -P 1    --nt=4
./prep.py -T 310 -P 1.05 --nt=4
...
./prep.py -T 410 -P 1.15 --nt=4
```
Here, `--nt=4` limits the number of threads to 4.

Locally, the above steps can be combined as
```
./mkTP.sh
```

For the second set, use `mkTP2.sh`.
For the two-dimensional set, use `mkTP2d.sh`.





### Run simulations ###

On lonestar
```
qsub TP_lonestar.pbs
```

Locally
```
cd T300P1
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx mdrun -v -deffnm npt
```



### Export energy ###

```
./TP_scan.sh
```
For the second set
```
./TP_scan.sh 2
```

Manually,
```
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx energy -f npt.edr -o ev.xvg
```
Select "Potential" (11), "Pressure" (16), and "Volume" (21).



### Run WHAM2 ###

```
cd ../../prog/gmx
ls --color=none -d ../../data/1VII/T[0-9][0-9][0-9]/ev.xvg > ev.ls
./xvgwham2 ev.ls --wham=mdiis
```

```
python xvgrun.py --ev
```



### Run MBAR2 ###

```
./xvgmbar2 ev.ls --wham=mdiis
```



### Benchmark ###

For the number of iteractions needed to reach an error tolerance
```
python xvgrun.py --ev
python stat.py
```

For the error versus the number of iteractions
```
python xvgtrace.py --ev
python stat.py
```





