# 0. Overview

This directory contains data for simulations on the villin headpiece.

Two figures in the manuscript depend on data collected here.

The first is the overall performance figure (Fig. 1).
This figure includes results from Ising model, Lennard-Jones fluid
and the villin headpiece, the system we have here.
This figure is drawn using the gnuplot script `nsnt.gp`.

The second figure is the error comparison figure (Fig. 5),
which is dedicated to the villin head piece.
This figure is drawn using the gnuplot script `whamcmp.gp`.

Below we shall describe how we prepared the system,
how simulation were run on the supercomputer,
and how the data were analysized.


1. System preparation
2. Running simulations
3. Generating energy log files
4. Running analyses
5. Making plots


# 1. System preparation

## 1.1 Master template under `init`

The initial system is saved in the directory `init`.
The files under that directory are the following

File        |  Description
------------|-------------------
1VII.pdb    | original PDB file
init.gro    | initial configuration
topol.top   | topology file
nvt.mdp     | template running parameters for NVT simulations
npt.mdp     | template running parameters for NPT simulations

In preparing the system, we went through the standard
GROMACS procedure.

The initial files were generated using the python script
`simulpdb.py`, which automates the GROMACS calls.


## 1.2 Using the template, `prep.py`

Once we have the template,
simulation files for different temperatures
were created using python `prep.py`.

We originally planned simulations for NVT and NPT ensembles,
however only data from the NVT ensemble simulations
were actually used.


## 1.3 NVT ensemble

The system under 12 temperatures, 300K, 310K, ... 410K
were prepared.

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

We have also prepared another set of simulations
for temperatures from 420K to 530K.
```
./mkT2.sh
```

## 1.4 NPT ensemble

We have created several sets of simulations
in the NPT ensemble.

To prepare systems one-by-one:
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




# 2. Running simulations

The simulations are run on the Lonestar supercomputer
To copy the necessary files

```
rsync -avz T[0-9][0-9][0-9] oo1@lonestar.tacc.utexas.edu:/scratch/02464/oo1/wham/data/1VII/
rsync -avz *.pbs            oo1@lonestar.tacc.utexas.edu:/scratch/02464/oo1/wham/data/1VII/
rsync -avz *.sh             oo1@lonestar.tacc.utexas.edu:/scratch/02464/oo1/wham/data/1VII/
```


## NVT ensemble simulations

On Lonestar
```
cds wham/data/1VII
module load python
qsub T_lonestar.pbs
```
Note, in order to submit jobs on Lonestar,
the module `python` must be loaded first.
The default version of `python` is too old.

Locally
```
cd T300
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx mdrun -v -deffnm nvt
```

## NPT ensemble

On lonestar
```
cds wham/data/1VII
module load python
qsub TP_lonestar.pbs
```

Locally
```
cd T300P1
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx mdrun -v -deffnm npt
```





# 3. Generating energy log files


## NVT ensemble

This is done by the shell script `T_scan.sh`
```
./T_scan.sh
```
Internally, it call the GROMACS program `gmx energy`.
For each directory T300, T310, ...,
it will create an energy log file `e.xvg` under it.

Preferrably, we should run this script on `lonestar`
so the `e.xvg` will be copied to the local machine
when we are running the syncin script.

If this process is done locally, then
```
./syncin && ./T_scan.sh
```

To manually run `gmx energy`,
do the following
```
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx energy -f nvt.edr -o e.xvg
```
Select "Potential" (11) and possibly "Pressure" (17).


## NPT ensemble

```
./TP_scan.sh
```
For the second set
```
./TP_scan.sh 2
```

Again if this process is done locally,
the first call
```
./syncin
```

Manually,
```
~/lwork/gmx/gromacs5.0/buildgcc/bin/gmx energy -f npt.edr -o ev.xvg
```
Select "Potential" (11), "Pressure" (16), and "Volume" (21).


## Generating lists of energy files

Data analysis programs such as `xvgwham` in `prog/gmx`,
requires a list of energy files as the input.
An example is `e.ls` here.
```
T300/e.xvg
T310/e.xvg
T320/e.xvg
T330/e.xvg
T340/e.xvg
T350/e.xvg
T360/e.xvg
T370/e.xvg
T380/e.xvg
T390/e.xvg
T400/e.xvg
T410/e.xvg
```

One way to generate these files is to use the `ls` command
```
ls --color=none -d T[0-9][0-9][0-9]/e.xvg > e.ls
```





# 4. Running analyses


Analyses can be done both locally and on Lonestar.

For Figure 1, the timing data are better run remotely.

For Figure 5, the error comparison data can be done locally.



## For Fig. 1, `nsnt.gp`


### Copy files

First copy the programs and analyses scripts to Lonestar.
```
make -C .. lonestar
```
Or manually,
```
rsync -avz ../*.pbs       oo1@lonestar.tacc.utexas.edu:/scratch/02464/oo1/wham/data/
rsync -avz ../dotdote.ls  oo1@lonestar.tacc.utexas.edu:/scratch/02464/oo1/wham/data/
```

### Prepare the energy log files

Run
```
cds wham/data/1VII
./T_scan.sh
```


### Running scripts

To submit the analysis job
```
cds wham/data
module load python
qsub 1VII.pbs
```


## For Fig. 5 `whamcmp.gp`

### Generate the reference values

The reference values of different methods are generated
using the following script
```
./mkwhamcmp.sh
```
In computing the references, all data points are used.
Particularlly, `mbar.out` is used as the reference
for other subsample data.


### Generate subsample data

```
cd whamcmpr0.01   && ./gen.sh && cd ..
cd whamcmpr0.0001 && ./gen.sh && cd ..
```

Note, if the data set changes, make sure to clear all log files
under `whamcmpr0.01` and `whamcmpr0.0001`



# 5. Making plots

Go to `doc/fig` (or `doc/figclr` for color figures).

## For Figure 1

```
gnuplot nsnt.gp
make nsnt.pdf
```

## For Figure 5

```
gnuplot whamcmp.gp
make whamcmp.pdf
```


# Files


File            | Description
----------------|--------------------------------------------------------------------
mkwhamcmp.sh    | script to generate files for the GNUplot script doc/fig/whamcmp.gp. This is sufficient only the old whacmp.gp, which applies to the entire data set. For the new version, which is compute for subsamples, check the subdirectories, whamcmpr0.01 and whamcmpr0.0001.  However, mbar.out is still needed, which can be generated either by mkwhamcmp.sh or manually.
xvgrun.py       | running script to measure the number of iterations and run time, output used for doc/fig/nsnt.gp, the output is .log
xvgtrace.py     | running script to measure the decay of error, not used, the output is .tr
xvgerr.py       | running script to estimate of error, output used for doc/fig/whamcmp.gp, the output is .log
stat.py         | script to analyze the data generated by xvgrun.py and xvgtrace.py
addlnzref.py    | add a column of the reference values of partition functions to whamxxx.out, used in `mkwhampcmp.sh`
T\_scan.sh      | generate e.xvg for NVT simulation trajectories
TP\_scan.sh     | generate ev.xvg for NPT simulation trajectories
mkT.sh          | prepare systems for NVT simulations
mkTPxxx.sh      | prepare systems for NPT simulations
prep.py         | prepare system for GROMACS simulations (low level), called in mkTxxx.sh and mkTPxxx.sh
T\_extend.sh    | extend simulation time (no longer needed)
uploadTPxxx.sh  | update initial files for running on Lonstar



# Usage of `xvgrun.py` and `xvgtrace.py`


`xvgrun.py` computes the number of iterations and run time.
`xvgtrace.py` gives a profile of the error versus the number of iterations.


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





