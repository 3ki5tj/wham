#!/usr/bin/env python



''' test WHAM for Lennard-Jones fluid '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 10000
nbases = 20
nequil = 5000
nsteps = 500000
fnlog = None
doev = False
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  WHAM on the Lennard-Jones fluid

  OPTIONS:
    -N              set the number of samples
    -M, --nbases=   set the maximal number of bases
    -m, --nequil=   set the number of equilibration steps
    -n, --nsteps=   set the number of simulation steps
    -o, --log=      set the output log file
    --ev, --lj2     do the two-dimensional case
    -v              be verbose
    --verbose=      set verbocity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvN:M:m:n:o:",
        [ "help", "verbose=",
          "nbases=", "nequil=", "nsteps=", "log=",
          "ev", "lj2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, nbases, nequil, nsteps, fnlog, doev, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-N",):
      nsamp = int(a)
    elif o in ("-M", "--nbases"):
      nbases = int(a)
    elif o in ("-m", "--nequil"):
      nequil = int(a)
    elif o in ("-n", "--nsteps"):
      nsteps = int(a)
    elif o in ("-o", "--log="):
      fnlog = a
    elif o in ("--ev", "--lj2"):
      doev = True
    elif o in ("-h", "--help"):
      usage()



def getnsteps(err):
  ''' get the nsteps from the error output '''

  ln = err.strip().split("\n")[-1]
  m = re.search(" ([0-9]+) steps", ln)
  if not m:
    print "Bad last line: ", err
    raise Exception
  return m.group(1) # keep it as a string, not an integer



def main():
  global fnlog

  zcom.runcmd("make -C ../../prog/lj")

  if doev:
    prog = "ljwham2"
    if not fnlog: fnlog = "lj2.log"
  else:
    prog = "ljwham"
    if not fnlog: fnlog = "lj.log"

  shutil.copy("../../prog/lj/" + prog, "./" + prog)

  cmd0 = "./%s --re --nequil=%d --nsteps=%d" % (
      prog, nequil, nsteps)

  ns = [0]*(nbases + 1)
  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    ret, out, err = zcom.runcmd(cmd0, capture = True)
    ns[0] = getnsteps(err)

    for nb in range(1, nbases + 1):
      cmd = "%s --wham=MDIIS --nbases=%d -H" % (cmd0, nb)
      ret, out, err = zcom.runcmd(cmd, capture = True)
      ns[nb] = getnsteps(err)

    # save to the log file
    open(fnlog, "a").write(" ".join(ns) + "\n")



if __name__ == "__main__":
  doargs()
  main()

