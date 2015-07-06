#!/bin/sh

opt="-r 0.01"

./xvgerr.py $opt --wham=mdiis
./xvgerr.py $opt --wham=mdiis --opt="--de=100" -o "xvgerrde100.log"
./xvgerr.py $opt --wham=st
./xvgerr.py $opt --wham=st    --opt="--de=100" -o "xvgerrde100.log"
./xvgerr.py $opt --wham=ui
./xvgerr.py $opt --wham=mdiis --mbar --opt="--tol=1e-8"
./xvgerr.py $opt --mbar --est --opt="--tol=1e-8"
