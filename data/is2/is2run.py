#!/usr/bin/env python



''' test WHAM for two-dimensional Ising model '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 100
nbases = 20
nequil = 100000
nsteps = 10000000
fnlog = "is2.log"
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  WHAM on the two-dimensional Ising model

  OPTIONS:
    -N              set the number of samples
    -M, --nbases=   set the maximal number of bases
    -m, --nequil=   set the number of equilibration steps
    -n, --nsteps=   set the number of simulation steps
    -o, --log=      set the output log file
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
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, nbases, nequil, nsteps, fnlog, verbose

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
    elif o in ("-h", "--help"):
      usage()



def getnsteps(err):
  ''' get the nsteps from the error output '''

  ln = err.strip().split("\n")[-1]
  m = re.search(" ([0-9]+) steps", ln)
  return m.group(1) # keep it as a string, not an integer



def main():
  zcom.runcmd("make -C ../../prog/is2")
  shutil.copy("../../prog/is2/is2wham", "./is2wham")

  cmd0 = "./is2wham --re --nequil=%d --nsteps=%d" % (
      nequil, nsteps)

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
