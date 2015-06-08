#!/usr/bin/env python



''' test WHAM for XVG files '''



import sys, os, getopt, shutil, re
import zcom



nsamp = 10000
radd = 0.1
nbases = 20
fnls = None
fnlog = None
update_method = " "
mthreshold = 10.0
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
    -r              set the sampling rate
    -M, --nbases=   set the maximal number of bases
    -l, --ls=       set the input list file
    -o, --log=      set the output log file
    --kth           use the KTH scheme in MDIIS
    --hp            use the HP scheme in MDIIS
    --mthreshold=   set the clean up threshold for MDIIS
    --opt=          set options to be passed to the command line
    --ev, --xvg2    do the two-dimensional case
    -v              be verbose
    --verbose=      set verbocity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvN:r:M:l:o:",
        [ "help", "verbose=",
          "nbases=", "KTH", "kth", "HP", "hp", "mthreshold=",
          "ls=", "log=",
          "opt=", "ev", "xvg2", "gmx2",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global nsamp, radd, nbases, update_method, mthreshold
  global fnls, fnlog, doev, cmdopt, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
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
    elif o in ("--mthreshold",):
      mthreshold = float(a)
    elif o in ("--opt",):
      cmdopt = a
    elif o in ("-l", "--ls="):
      fnls = a
    elif o in ("-o", "--log="):
      fnlog = a
    elif o in ("--ev", "--xvg2", "--gmx2"):
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
  global radd, cmdopt, fnls, fnlog

  zcom.runcmd("make -C ../../prog/gmx")

  if doev:
    prog = "xvgwham2"
    if not fnlog: fnlog = "xvg2.log"
    if not fnls: fnls = "ev.ls"
  else:
    prog = "xvgwham"
    if not fnlog: fnlog = "xvg.log"
    if not fnls: fnls = "e.ls"

  shutil.copy("../../prog/gmx/" + prog, "./" + prog)

  cmd0 = "./%s -r %g %s %s" % (
      prog, radd, cmdopt, fnls)

  ns = [0]*(nbases + 1)
  for i in range(nsamp):
    print "running sample %d/%d..." % (i, nsamp)

    # use the direct WHAM
    ret, out, err = zcom.runcmd(cmd0, capture = True)
    ns[0] = getnsteps(err)

    for nb in range(1, nbases + 1):
      cmd = "%s --wham=MDIIS --nbases=%d -H %s --mthreshold=%g" % (
          cmd0, nb, update_method, mthreshold)
      ret, out, err = zcom.runcmd(cmd, capture = True)
      ns[nb] = getnsteps(err)

    # save to the log file
    open(fnlog, "a").write(" ".join(ns) + "\n")



if __name__ == "__main__":
  doargs()
  main()

