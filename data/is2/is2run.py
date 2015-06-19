#!/usr/bin/env python



''' test WHAM for two-dimensional Ising model
    find the number of iterations needed to reach an error tolerance '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 10000
nbases = 20
nequil = " "
nsteps = " "
fnlog = "is2.log"
update_method = " "
mthreshold = " "
tol = " "
cmdopt = ""
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  WHAM on the two-dimensional Ising model

  OPTIONS:
    -N              set the number of samples
    -M, --nbases=   set the maximal number of bases in MDIIS
    -m, --nequil=   set the number of equilibration steps
    -n, --nsteps=   set the number of simulation steps
    -o, --log=      set the output log file
    --kth           use the KTH scheme in MDIIS
    --hp            use the HP scheme in MDIIS
    --hpl           use the HPL scheme in MDIIS
    --mthreshold=   set the cleanup threshold for MDIIS
    --tol=          set the tolerance of error
    --opt=          set options to be passed to the command line
    -v              be verbose
    --verbose=      set verbosity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvN:M:m:n:o:",
        [ "help", "verbose=",
          "nbases=", "KTH", "kth", "HP", "hp", "HPL", "hpl",
          "mthreshold=", "tol=",
          "nequil=", "nsteps=", "log=", "opt=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, nbases, update_method, mthreshold, tol
  global nequil, nsteps, fnlog, cmdopt, verbose

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
    elif o in ("--HPL", "--hpl"):
      update_method = "--hpl"
    elif o in ("--mthreshold",):
      mthreshold = "--mthreshold=%g" % float(a)
    elif o in ("--tol",):
      tol = "--tol=%g" % float(a)
    elif o in ("--opt",):
      cmdopt = a
    elif o in ("-m", "--nequil"):
      nequil = "--nequil=%d" % int(a)
    elif o in ("-n", "--nsteps"):
      nsteps = "--nsteps=%d" % int(a)
    elif o in ("-o", "--log="):
      fnlog = a
    elif o in ("-h", "--help"):
      usage()



def getnstepstime(err):
  ''' get the nsteps and time from the error output '''

  ln = err.strip().split("\n")[-1].strip()

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
  fnhis = arr[0] + "_hist.dat"

  progdir = "../../prog"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir

  zcom.runcmd("make -C %s/is2" % progdir)

  prog = "is2wham"

  shutil.copy("%s/is2/%s" % (progdir, prog), "./%s" % prog)

  cmd0 = "./%s --fnhis=%s %s %s" % (
      prog, fnhis, tol, cmdopt)
  cmd0 = cmd0.strip()

  ns = [0]*(nbases + 1)
  tm = [0]*(nbases + 1)
  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    cmd = "%s --re %s %s" % (cmd0, nequil, nsteps)
    ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
    ns[0], tm[0] = getnstepstime(err)

    for nb in range(1, nbases + 1):
      cmd = "%s --wham=MDIIS --nbases=%d -H %s %s" % (
          cmd0, nb, update_method, mthreshold)
      ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
      ns[nb], tm[nb] = getnstepstime(err)

    # save to the log files
    open(fnlog, "a").write(" ".join(ns) + "\n")
    open(fntmlog, "a").write(" ".join(tm) + "\n")



if __name__ == "__main__":
  doargs()
  main()

