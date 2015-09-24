#!/usr/bin/env python



''' test WHAM for XVG files
    find the number of iterations need to reach an error tolerance '''



import sys, os, getopt, shutil, re
import zcom



mbar = False
nsamp = 10000
radd = 1.0
nbases = 20
fnls = None
fnlog = None
fnact = "--fnact=act.dat"
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
  WHAM on XVG files

  OPTIONS:
    --mbar          set the program to MBAR
    -N              set the number of samples
    -r              set the sampling rate
    -M, --nbases=   set the maximal number of bases
    -l, --ls=       set the input list file
    -o, --log=      set the output log file
    --kth           use the KTH scheme in MDIIS
    --hp            use the HP scheme in MDIIS
    --hpl           use the HP scheme in MDIIS
    --mthreshold=   set the clean up threshold for MDIIS
    --tol=          set the tolerance of error
    --opt=          set options to be passed to the command line
    --ev, --xvg2    do the two-dimensional case
    -v              be verbose
    --verbose=      set verbosity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvN:r:M:l:o:",
        [ "help", "verbose=",
          "mbar",
          "nbases=", "KTH", "kth", "HP", "hp", "HPL", "hpl",
          "mthreshold=", "tol=",
          "ls=", "log=",
          "opt=", "ev", "xvg2", "gmx2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global mbar, nsamp, radd, nbases, update_method, mthreshold, tol
  global fnls, fnlog, doev, cmdopt, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("--mbar",):
      mbar = True
    elif o in ("-N",):
      nsamp = int(a)
    elif o in ("-r",):
      radd = float(a)
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
    elif o in ("-l", "--ls"):
      fnls = a
    elif o in ("-o", "--log"):
      fnlog = a
    elif o in ("--ev", "--xvg2", "--gmx2"):
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
  global mbar, radd, cmdopt, fnls, fnlog, fnact

  progdir = "../../prog"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir

  zcom.runcmd("make -C %s/gmx" % progdir)

  if doev:
    if mbar:
      prog = "xvgmbar2"
    else:
      prog = "xvgwham2"
    if not fnlog: fnlog = "xvg2.log"
    if not fnls: fnls = "ev.ls"
    hisopt = "--fnhis2"
    fnhis = "hist2.dat"
  else:
    if mbar:
      prog = "xvgmbar"
    else:
      prog = "xvgwham"
    if not fnlog: fnlog = "xvg.log"
    if not fnls: fnls = "e.ls"
    hisopt = "--fnhis"
    fnhis = "hist.dat"

  if mbar:
    strfnhis = ""
    # modify the the log file name
    arr = os.path.splitext(fnlog)
    fnlog = arr[0] + "_mbar" + arr[1]
  else:
    arr = os.path.splitext(fnlog)
    fnlog = arr[0] + "_wham" + arr[1]
    fnhis = arr[0] + "_" + fnhis
    strfnhis = hisopt + "=" + fnhis

  arr = os.path.splitext(fnlog)
  fntmlog = arr[0] + "tm" + arr[1]

  try:
    shutil.copy("%s/gmx/%s" % (progdir, prog), "./%s" % prog)
  except:
    pass

  cmd0 = "./%s %s %s %s %s %s" % (
      prog, strfnhis, fnls, tol, fnact, cmdopt)
  cmd0 = cmd0.strip()

  if radd < 1:
    cmd0 += " -r %g" % radd

  if not mbar:
    # run for the first time to save the histogram
    zcom.runcmd(cmd0 + " --wham=ST")

    # we will load the histogram, and bootstrap later on
    cmd0 += " --bootstrap -H"
  else:
    # no need to bootstrap if only a random subsample is used
    if radd >= 1:
      cmd0 += " --bootstrap"

  cmd0 = cmd0.strip()

  ns = [0]*(nbases + 1)
  tm = [0]*(nbases + 1)
  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    ret, out, err = zcom.runcmd(cmd0, capture = True)
    ns[0], tm[0] = getnstepstime(err)

    for nb in range(1, nbases + 1):
      cmd = "%s --%s=MDIIS --nbases=%d %s %s" % (
          cmd0, "mbar" if mbar else "wham",
          nb, update_method, mthreshold)
      ret, out, err = zcom.runcmd(cmd.strip(), capture = True)
      ns[nb], tm[nb] = getnstepstime(err)

    # save to the log files
    open(fnlog, "a").write(" ".join(ns) + "\n")
    open(fntmlog, "a").write(" ".join(tm) + "\n")



if __name__ == "__main__":
  doargs()
  main()

