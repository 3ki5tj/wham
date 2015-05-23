## Run it ##

1. Generate `is2.log`
```
./is2run.py
```
This command can be run multiple times to increase the sample size.


2. Generate `is2wham.dat`
```
./stat.py
```

3. Plot the result in gnuplot
```
set logscale y
plot "is2wham.dat" u 1:2:3 w e pt 7
```

