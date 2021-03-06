#!/usr/bin/env python



''' add reference lnz from MBAR '''



import sys, os, getopt, shutil, re, glob


fninp = None
fnref = "mbar.out"
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] file""" % sys.argv[0]

  print """
  OPTIONS:
    --ref=          set the reference file
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
        [ "ref=",
          "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fninp, fnref

  for o, a in opts:
    if o in ("--ref",):
      fnref = a
    elif o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help"):
      usage()

  if len(args) > 0:
    fninp = args
  else:
    fninp = glob.glob("*wham*.out") + glob.glob("*mbar*.out") + glob.glob("est*.out")



def main(fn, fnref):
  # load the reference array
  s = open(fnref).readlines()
  n = len(s)
  betref = [0]*n
  arrref = [0]*n
  for i in range(n):
    ln = s[i].strip()
    x = ln.split()
    try:
      betref[i] = float(x[1])
      # do not convert to float
      arrref[i] = x[2]
    except ValueError:
      print ln, x
      raw_input()

  # load the input data
  s = open(fn).readlines()
  nn = len( [ln for ln in s if not ln.startswith("#")] )
  if nn != n:
    print "number of lines mismatch %s(%s) vs %s " % (nn, fn, n)
    raise Exception
  bet = [-0]*n
  arr = [-0]*n

  maxcol = 6
  if fn.startswith("est"):
    maxcol = 11
  ii = 0
  for i in range(n):
    # find the ith noncomment line
    while s[ii].startswith("#"):
      ii += 1
    ln = s[ii].rstrip()
    arr = ln.split()
    if len(arr) > maxcol:
      # assuming every column beyond column `maxcol` is added
      # by addlnzref, so we can remove them safely
      p = ln.rfind( arr[maxcol] )
      if p >= 0:
        ln = ln[:p].rstrip()
    s[ii] = ln + "\t" + arrref[i] + "\n"
    ii += 1

  print "updating %s" % fn
  open(fn, "w").writelines(s)



if __name__ == "__main__":
  doargs()
  for fn in fninp:
    if fn != fnref:
      main(fn, fnref)

