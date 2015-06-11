
## Reference ##

`histref.dat` was used to produce
`is2nb0.out` and `is2nb10.out`.


## Number of iterations versus number of bases ##

1. Generate `is2.log` and `is2tm.log`
```
./is2run.py
```
This command can be run multiple times to increase the sample size.


2. Generate `is2wham.dat` (for the number of steps)
and `is2tmwham.dat` (for the runtime)
```
./stat.py
```

3. Plot the result in gnuplot
```
set logscale y
plot "is2wham.dat" u 1:2:3 w e pt 7
plot "is2tmwham.dat" u 1:2:3 w e pt 7
```


## Error versus number of iterations ##

Use `is2trace.py`
