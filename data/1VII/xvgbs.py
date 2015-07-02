#!/usr/bin/env python



''' test WHAM for XVG files
    estimate the error  '''



import sys, os, getopt, shutil, re
import zcom



mbar = False
nsamp = 10000
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
    --mbar          set the program to MBAR
    -N              set the number of samples
    -M, --nbases=   set the number of bases
    -l, --ls=       set the input list file
    -o, --log=      set the output log file
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
        "hvN:M:l:o:",
        [ "help", "verbose=",
          "mbar",
          "nbases=",
          "ls=", "log=",
          "opt=", "ev", "xvg2", "gmx2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global mbar, nsamp
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
  global mbar, cmdopt, fnls, fnlog

  progdir = "../../prog"
  if not os.path.isdir(progdir):
    progdir = "../" + progdir

  zcom.runcmd("make -C %s/gmx" % progdir)

  if doev:
    if mbar:
      prog = "xvgmbar2"
    else:
      prog = "xvgwham2"
    if not fnlog: fnlog = "xvg2bs.log"
    if not fnls: fnls = "ev.ls"
    hisopt = "--fnhis2"
    fnhis = "histbs2.dat"
  else:
    if mbar:
      prog = "xvgmbar"
    else:
      prog = "xvgwham"
    if not fnlog: fnlog = "xvgbs.log"
    if not fnls: fnls = "e.ls"
    hisopt = "--fnhis"
    fnhis = "histbs.dat"

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

  try:
    shutil.copy("%s/gmx/%s" % (progdir, prog), "./%s" % prog)
  except:
    pass

  cmd0 = "./%s -v %s %s %s --%s=MDIIS %s" % (
      prog, strfnhis, fnls, cmdopt,
      "mbar" if mbar else "wham", nbases)
  cmd0 = cmd0.strip()

  if not mbar:
    # run for the first time to save the histogram
    zcom.runcmd(cmd0)
    # we will load the histogram, and bootstrap later on
    cmd0 += " --bootstrap -H"
  else:
    cmd0 += " --bootstrap"

  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    ret, out, err = zcom.runcmd(cmd0.strip(), capture = True)
    if i == 0:
      bet = getcol(out, 1)
      open(fnlog, "a").write("#HEADER " + " ".join(bet) + "\n")
    lnz = getcol(out, 2)

    # save to the log files
    open(fnlog, "a").write(" ".join(lnz) + "\n")



if __name__ == "__main__":
  doargs()
  main()

