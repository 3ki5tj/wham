#!/usr/bin/env python



import os, sys, getopt, shutil, re, glob
from math import *



geomean = -1
epsilon = 1e-10
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] input.log""" % sys.argv[0]

  print """
  Compute statistics

  OPTIONS:
    --geo           geometric mean
    --eps=          set the minimal value
    -v              be verbose
    --verbose=      set verbocity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hv",
        [ "geo", "eps=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global geomean, epsilon, verbose

  for o, a in opts:
    if o in ("--geo",):
      geomean = 1
    elif o in ("--eps",):
      epsilon = float(a)
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("-h", "--help"):
      usage()

  if len(args):
    return args
  else:
    return glob.glob("*.log") + glob.glob("*.tr")



def dostat(fninp):
  global geomean, epsilon

  if geomean < 0: # default
    arr = os.path.splitext(fninp)
    if arr[1] == ".tr":
      geom = True
    else:
      geom = False
  else:
    geom = geomean

  s = open(fninp).readlines()

  n = len(s)

  # find the longest length
  m = 0
  for i in range(n):
    m = max(m, len(s[i].strip().split()))

  sy = [0.0] * m
  syy = [0.0] * m
  for i in range(n):
    # parse the line to an array of floating-point numbers
    arr = [float(x) for x in s[i].strip().split()]
    for j in range(m):
      if j >= len(arr):
        if geom:
          x = epsilon
        else:
          x = 0
      else:
        x = arr[j]
      if geom:
        x = log(x)
      sy[j] += x
      syy[j] += x * x

  ave = [0.0] * m
  var = [0.0] * m
  for j in range(m):
    ave[j] = sy[j] / n;
    var[j] = syy[j] / n - ave[j] * ave[j]
    if geom:
      ave[j] = exp( ave[j] )

  # write the output file
  s = "# %d %d\n" % (n, m)
  for j in range(m):
    if n >= 2:
      std = (var[j]*n/(n-1)) ** 0.5
    else:
      std = 0
    s += "%d %g %g\n" % (j, ave[j], std)

  fnout = fninp.split(".")[0] + "wham.dat"
  open(fnout, "w").write(s)
  print "save results to %s" % fnout



if __name__ == "__main__":
  fninps = doargs()
  for fninp in fninps:
    dostat(fninp)
