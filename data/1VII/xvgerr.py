#!/usr/bin/env python



''' test WHAM for XVG files
    estimate the error '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 10000
bootstrap = ""
radd = "-r 0.0001"
est = False
mbar = False
whammethod = "MDIIS"
fnls = None
fnlog = None
nbases = " "
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
    -N              set the number of samples
    --est           compute rough estimations, e.g., average & BAR
    --mbar          set the program to MBAR
    -r, --radd=     set the rate of adding frames
    --bootstrap     use bootstrapping
    --wham=         WHAM or MBAR method, MDIIS, ST, UI
    -M, --nbases=   set the number of bases in MDIIS
    -l, --ls=       set the input list file
    -o, --log=      set the output log file, e.g., xvgerr.log
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
          "radd=", "bootstrap", "est", "mbar",
          "wham=", "nbases=",
          "ls=", "log=",
          "opt=", "ev", "xvg2", "gmx2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, radd, bootstrap, est, mbar, whammethod
  global fnls, fnlog, doev, cmdopt, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-N",):
      nsamp = int(a)
    elif o in ("-r", "--radd"):
      radd = "-r %g" % float(a)
    elif o in ("--bootstrap"):
      bootstrap = "--bootstrap"
    elif o in ("--est",):
      est = True
      mbar = True  # estimations are done through the mbar program
    elif o in ("--mbar",):
      mbar = True
    elif o in ("--wham",):
      whammethod = a
    elif o in ("-M", "--nbases"):
      nbases = "--nbases=%d" % int(a)
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




def getcol(out, col):
  ''' get the partition from the output '''

  s = out.strip().split("\n")
  n = len(s)
  arr = [""]*n

  for i in range(n):
    ln = s[i].strip()
    arr[i] = ln.split()[col]

  return arr



def main():
  global bootstrap, est, mbar, cmdopt, fnls, fnlog

  progdir = "../../prog"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir

  zcom.runcmd("make -C %s/gmx" % progdir)

  if doev:
    if mbar:
      prog = "xvgmbar2"
    else:
      prog = "xvgwham2"
    if not fnlog: fnlog = "xvg2err.log"
    if not fnls: fnls = "ev.ls"
    hisopt = "--fnhis2"
    fnhis = "histerr2.dat"
  else:
    if mbar:
      prog = "xvgmbar"
    else:
      prog = "xvgwham"
    if not fnlog: fnlog = "xvgerr.log"
    if not fnls: fnls = "e.ls"
    hisopt = "--fnhis"
    fnhis = "histerr.dat"

  if mbar:
    strfnhis = ""
    # modify the the log file name
    arr = os.path.splitext(fnlog)
    if not est:
      fnlog = arr[0] + "_mbar" + arr[1]
  else:
    arr = os.path.splitext(fnlog)
    wmethod = whammethod.upper()
    if wmethod.startswith("UI"):
      method = "ui"
    elif wmethod == "ST":
      method = "st"
    else:
      method = ""
    fnlog = arr[0] + "_wham" + method + arr[1]
    fnhis = arr[0] + method + "_" + fnhis
    strfnhis = hisopt + "=" + fnhis

  try:
    shutil.copy("%s/gmx/%s" % (progdir, prog), "./%s" % prog)
  except:
    pass

  cmd0 = "./%s -v %s %s %s %s %s %s --%s=%s %s" % (
      prog, strfnhis, fnls, cmdopt, radd, bootstrap,
      "--est" if est else "",
      "mbar" if mbar else "wham", whammethod, nbases)
  cmd0 = cmd0.strip()

  if bootstrap and not mbar:
    # run for the first time to save the histogram
    zcom.runcmd(cmd0)
    # we will load the histogram, and bootstrap later on
    cmd0 += " -H"

  if est:
    arr = os.path.splitext(fnlog)
    fnlogs = [
      arr[0] + "avea_mbar" + arr[1],
      arr[0] + "aveb_mbar" + arr[1],
      arr[0] + "avec_mbar" + arr[1],
      arr[0] + "expa_mbar" + arr[1],
      arr[0] + "expb_mbar" + arr[1],
      arr[0] + "bar_mbar" + arr[1],
      arr[0] + "gp_mbar" + arr[1],
      arr[0] + "tg_mbar" + arr[1],
      arr[0] + "lnv_mbar" + arr[1],
      ]
  else:
    fnlogs = [ fnlog ]

  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # run WHAM
    ret, out, err = zcom.runcmd(cmd0.strip(), capture = True)
    if i == 0:
      bet = getcol(out, 1)
      hdr = "#HEADER " + " ".join(bet) + "\n"
      for fnlog in fnlogs:
        open(fnlog, "a").write(hdr)

    for i in range(len(fnlogs)):
      fnlog = fnlogs[i]
      col = i + 2
      lnz = getcol(out, col)
      open(fnlog, "a").write(" ".join(lnz) + "\n")



if __name__ == "__main__":
  doargs()
  main()

