#!/usr/bin/env python



''' test WHAM for Lennard-Jones fluid '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 10000
nbases = 20
nequil = 5000
nsteps = 500000
fnlog = None
update_method = " "
mthreshold = " "
tol = " "
cmdopt = ""
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
    --kth           use the KTH scheme in MDIIS
    --hp            use the HP scheme in MDIIS
    --mthreshold=   set the clean up threshold for MDIIS
    --tol=          set the tolerance of error
    --opt=          set options to be passed to the command line
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
          "nbases=", "KTH", "kth", "HP", "hp",
          "mthreshold=", "tol=",
          "nequil=", "nsteps=", "log=",
          "opt=", "ev", "lj2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, nbases, update_method, mthreshold
  global nequil, nsteps, fnlog, doev, cmdopt, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-N",):
      nsamp = int(a)
    elif o in ("-M", "--nbases"):
      nbases = int(a)
    elif o in ("--KTH", "--kth"):
      update_method = "--kth"
    elif o in ("--HP", "--hp"):
      update_method = "--hp"
    elif o in ("--mthreshold",):
      mthreshold = "--mthreshold=%g" % float(a)
    elif o in ("--tol",):
      tol = "--tol=%g" % float(a)
    elif o in ("--opt",):
      cmdopt = a
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



def getnstepstime(err):
  ''' get the nsteps and time from the error output '''

  ln = err.strip().split("\n")[-1]

  m = re.search(" ([0-9]+) steps", ln)
  if not m:
    print "line %s: no number of iterations" % ln
    raise Exception
  ns = m.group(1) # keep it as a string, not an integer

  m = re.search("time ([0-9.]+)s", ln)
  if not m:
    print "line %s: no timing information" % ln
    raise Exception
  tm = m.group(1)

  return ns, tm



def main():
  global cmdopt, fnlog

  arr = os.path.splitext(fnlog)
  fntmlog = arr[0] + "tm" + arr[1]

  zcom.runcmd("make -C ../../prog/lj")

  if doev:
    prog = "ljwham2"
    if not fnlog: fnlog = "lj2.log"
  else:
    prog = "ljwham"
    if not fnlog: fnlog = "lj.log"

  shutil.copy("../../prog/lj/" + prog, "./" + prog)

  cmd0 = "./%s --re --nequil=%d --nsteps=%d %s %s" % (
      prog, nequil, nsteps, tol, cmdopt)

  ns = [0]*(nbases + 1)
  tm = [0]*(nbases + 1)
  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    ret, out, err = zcom.runcmd(cmd0, capture = True)
    ns[0], tm[0] = getnstepstime(err)

    for nb in range(1, nbases + 1):
      cmd = "%s --wham=MDIIS --nbases=%d -H %s %s" % (
          cmd0.strip(), nb, update_method, mthreshold)
      ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
      ns[nb], tm[nb] = getnstepstime(err)

    # save to the log files
    open(fnlog, "a").write(" ".join(ns) + "\n")
    open(fntmlog, "a").write(" ".join(tm) + "\n")



if __name__ == "__main__":
  doargs()
  main()

