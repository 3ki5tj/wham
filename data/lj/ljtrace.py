#!/usr/bin/env python



''' test WHAM for Lennard-Jones fluid
    trace the error versus the number of iterations '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 10000
nbases = 20
dnbases = 5
nequil = 5000
nsteps = 500000
fntr = None
update_method = " "
mthreshold = " "
itmin = "--itmin=100"
itmax = "--itmax=1000"
tol = "--tol=1e-10"
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
    -M, --nbases=   set the maximal number of bases in MDIIS
    -D, --dnbases=  set the step size of the number of bases in MDIIS
    -m, --nequil=   set the number of equilibration steps
    -n, --nsteps=   set the number of simulation steps
    -o, --trace=    set the output trace file
    --kth           use the KTH scheme in MDIIS
    --hp            use the HP scheme in MDIIS
    --hpl           use the HPL scheme in MDIIS
    --mthreshold=   set the clean up threshold for MDIIS
    --itmin=        set the minimal number of iterations
    --itmax=        set the maximal number of iterations
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
        "hvN:M:D:m:n:o:",
        [ "help", "verbose=",
          "nbases=", "dnbases=", "KTH", "kth", "HP", "hp", "HPL", "hpl",
          "mthreshold=", "itmin=", "itmax=", "tol=",
          "nequil=", "nsteps=", "trace=",
          "opt=", "ev", "lj2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, nbases, dnbases, update_method
  global mthreshold, itmin, itmax, tol
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
    elif o in ("-D", "--dnbases"):
      dnbases = int(a)
    elif o in ("--KTH", "--kth"):
      update_method = "--kth"
    elif o in ("--HP", "--hp"):
      update_method = "--hp"
    elif o in ("--HPL", "--hpl"):
      update_method = "--hpl"
    elif o in ("--mthreshold",):
      mthreshold = "--mthreshold=%g" % float(a)
    elif o in ("--itmin",):
      itmin = "--itmin=%d" % int(a)
    elif o in ("--itmax",):
      itmax = "--itmax=%d" % int(a)
    elif o in ("--tol",):
      tol = "--tol=%g" % float(a)
    elif o in ("--opt",):
      cmdopt = a
    elif o in ("-m", "--nequil"):
      nequil = int(a)
    elif o in ("-n", "--nsteps"):
      nsteps = int(a)
    elif o in ("-o", "--trace="):
      fntr = a
    elif o in ("--ev", "--lj2"):
      doev = True
    elif o in ("-h", "--help"):
      usage()



def gettrace(nb, err):
  ''' get the errors vs. step from the output '''

  s = err.strip().split("\n")

  tr = {}
  for ln in s:
    if not ln.startswith("it "):
      continue

    m = re.search("it ([0-9]+),", ln)
    if not m:
      print "line %s: no number of iterations" % ln
      raise Exception
    it = int( m.group(1) )

    m = re.search("err[a-z]* ([0-9.e+-]+) -> ([0-9.e+-]+)", ln)
    if not m:
      print "line %s: no error information" % ln
      raise Exception
    err0 = float( m.group(1) )
    err = float( m.group(2) )

    if it == 1:
      tr[0] = err0
    tr[it] = err

  itmax = max(k for k in tr)
  s = ""
  for k in range(itmax + 1):
    s += "%g " % tr[k]
  s = s.strip()

  arr = os.path.splitext(fntr)
  fntrnb = "%s_nb%s%s" % (arr[0], nb, arr[1])
  open(fntrnb, "a").write(s + "\n")
  return tr



def main():
  global cmdopt, fntr

  zcom.runcmd("make -C ../../prog/lj")

  if doev:
    prog = "ljwham2"
    if not fntr: fntr = "lj2.tr"
    fnhis = "hist2.dat"
  else:
    prog = "ljwham"
    if not fntr: fntr = "lj.tr"
    fnhis = "hist.dat"

  arr = os.path.splitext(fntr)
  fnhis = arr[0] + "_tr_" + fnhis

  try:
    shutil.copy("../../prog/lj/%s" % prog, "./%s" % prog)
  except:
    pass

  cmd0 = "./%s -v %s %s %s %s" % (
      prog, itmin, itmax, tol, cmdopt)
  cmd0 = cmd0.strip()

  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    cmd = "%s --re --nequil=%d --nsteps=%d" % (cmd0, nequil, nsteps)
    ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
    gettrace(0, err)

    for nb in range(dnbases, nbases + 1, dnbases):
      cmd = "%s --wham=MDIIS --nbases=%d -H %s %s" % (
          cmd0, nb, update_method, mthreshold)
      ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
      gettrace(nb, err)



if __name__ == "__main__":
  doargs()
  main()

