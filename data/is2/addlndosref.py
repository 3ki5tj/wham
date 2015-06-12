#!/usr/bin/env python



''' add reference density of states '''



import sys, os, getopt, shutil, re, glob


fninp = None
fnref = "islogdos64x64.txt"
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS] file""" % sys.argv[0]

  print """
  OPTIONS:
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
        [ "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fninp, fnref

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-h", "--help"):
      usage()

  if len(args) > 0:
    fninp = args
  else:
    fninp = glob.glob("lndos*.dat")



def main(fn, fnref):
  # load the reference array
  s = open(fnref).read().strip().split('\n')
  n = len(s) - 1
  arrref = [0]*(n + 1)
  for i in range(n + 1):
    ln = s[i]
    x = ln.strip()
    p = x.find("`")
    if p >= 0: x = x[:p]
    try:
      arrref[i] = float(x)
    except ValueError:
      print ln, x
      raw_input()

  # load the input data
  s = open(fn).readlines()
  arr = [-10000]*(n + 1)
  imin = 100000000
  imax = 0
  for ln in s:
    xy = ln.strip().split()
    x = int(xy[0])
    y = float(xy[1])
    i = (x + 2*n) / 4
    imin = min(i, imin)
    imax = max(i, imax)
    arr[i] = y

  # align arr[] with arrref[]
  dy = 0
  cnt = 0
  for i in range(imin, imax + 1):
    if arrref[i] < -100 or arr[i] < -100:
      continue
    dy += arrref[i] - arr[i]
    cnt += 1
  dy /= cnt

  s = ""
  for i in range(imin, imax + 1):
    arr[i] += dy
    s += "%d %s %s\n" % (-2*n+i*4, arr[i], arrref[i])

  print "updating %s" % fn
  open(fn, "w").write(s)



if __name__ == "__main__":
  doargs()
  for fn in fninp:
    main(fn, fnref)

