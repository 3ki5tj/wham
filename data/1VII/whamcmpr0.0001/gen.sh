#!/bin/sh


./xvgerr.py --wham=st
./xvgerr.py --wham=st    --opt="--de=100" -o "xvgerrde100.log"
./xvgerr.py --wham=ui
./xvgerr.py --wham=ui    --opt="--de=100" -o "xvgerrde100.log"
./xvgerr.py --wham=mdiis
./xvgerr.py --wham=mdiis --opt="--de=100" -o "xvgerrde100.log"
./xvgerr.py --mbar --est --opt="--tol=1e-8"
./xvgerr.py --wham=mdiis --mbar --opt="--tol=1e-8"

make err
