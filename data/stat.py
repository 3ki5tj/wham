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
        [ "geo=", "eps=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global geomean, epsilon, verbose

  for o, a in opts:
    if o in ("--geo",):
      geomean = int(a)
    elif o in ("--eps",):
      epsilon = float(a)
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("-h", "--help"):
      usage()

  if len(args):
    return args
  else:
    return glob.glob("*.log") + glob.glob("*.tr") + glob.glob("*.bs")



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

  # number of lines, which is the sample size
  n = 0
  for i in range(len(s)):
    if not s[i].startswith("#"):
      n += 1

  s0 = s
  s = []
  header = None
  for i in range(len(s0)):
    if not s0[i].startswith("#"):
      s += [ s0[i], ]
    elif s0[i].startswith("#HEADER"):
      header = s0[i][7:].strip().split()

  # find the maximal number of columns
  m = 0
  for i in range(n):
    m = max(m, len(s[i].strip().split()))

  # compute the statistical moments for each column
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
    # compute the average and variances for column j
    ave[j] = sy[j] / n;
    var[j] = syy[j] - ave[j] * ave[j] * n
    if n > 1:
      var[j] /= n - 1
    if geom:
      ave[j] = exp( ave[j] )

  # write the output file
  s = "# %d %d\n" % (n, m)
  for j in range(m):
    hdr = header[j] if header else ""
    s += "%d %s %g %g\n" % (j, hdr, ave[j], var[j] ** 0.5)

  # construct the output file name
  # the input file name is usually xxx.log, xxx.tr, or xxx.bs
  nm = os.path.splitext(fninp)[0]
  p = nm.find("_wham")
  q = nm.find("_mbar")
  if p >= 0:
    fnout = nm[:p] + nm[p+5:] + "wham.dat"
  elif q >= 0:
    fnout = nm[:q] + nm[q+5:] + "mbar.dat"
  else:
    fnout = nm + "wham.dat"

  open(fnout, "w").write(s)
  print "save results to %s" % fnout



if __name__ == "__main__":
  fninps = doargs()
  for fninp in fninps:
    dostat(fninp)
