#!/usr/bin/env python




''' prepare a system for MD simulation
    for NVT and NPT ensembles
    for precise setting
'''


import sys, os, subprocess, getopt, shutil, re, random, glob, math
import zcom




temp = 300
pres = None
nequil = 10000
nsteps = 100000000
nthreads = 0

verbose = 0
srcroot = "$HOME/work/gmx/gromacs5.0"
buildroot = srcroot + "/buildiccdbl"
gmxtopdir = srcroot + "/share/top"



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  prepare a system for GROMACS MD simulation from a PDB file

  OPTIONS:
    -T              set the temperature
    -P              set the pressure
    -m              set the number of equilibration steps
    -n              set the number of steps
    --nt=           set the number of threads
    -v              be verbose
    --verbose=      set verbocity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "hvt:T:p:P:m:n:",
        [ "help", "verbose=",
          "nt=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global temp, pres, nequil, nsteps, nthreads, verbose
  global srcroot, buildroot, gmxtopdir

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-t", "-T",):
      temp = float(a)
    elif o in ("-p", "-P",):
      pres = float(a)
    elif o in ("-m",):
      nequil = int(a)
    elif o in ("-n",):
      nsteps = int(a)
    elif o in ("--nt",):
      nthreads = int(a)
    elif o in ("-h", "--help"):
      usage()

  ret, out, err = zcom.runcmd("uname -a", capture = True)
  if out.find("host.utmb.edu") >= 0:
    srcroot = "$HOME/lwork/gmx/gromacs5.0"
    buildroot = srcroot + "/buildgccdbl"
    gmxtopdir = srcroot + "/share/top"



def changensteps(fn, nsteps):
  ''' set the number of steps in the .mdp file '''
  s = open(fn).readlines()
  for i in range(len(s)):
    if s[i].startswith("nsteps"):
      s[i] = "nsteps = %d\n" % nsteps
  open(fn, "w").writelines(s)



def mkmdp(fn, T = None, P = None):
  ''' make .mdp file '''
  s = open(fn).readlines()
  if T: # set the temperature
    for i in range(len(s)):
      if s[i].startswith("ref_t") or s[i].startswith("gen_temp"):
        s[i] = s[i].replace("300", "%g" % T)

  if P: # set the pressure
    for i in range(len(s)):
      if s[i].startswith("ref_p"):
        s[i] = s[i].replace("1.0", "%g" % P)

  # change the random number seed
  for i in range(len(s)):
    if s[i].startswith("gen_seed"):
      s[i] = "gen_seed = %d\n" % random.randint(1, 2000000000)

  return s



def mklonestarpbs(simulname):
  ''' make foo.pbs for room temperature simulation
      this script is provided for convenience '''

  global srcroot, buildroot, gmxtopdir
  mdrun = buildroot + "/bin/gmx_d mdrun"

  d = {
      "simulname" : simulname,
      "mdrun" : mdrundir,
      "gmxtopdir" : gmxtopdir,
      }
  return zcom.templrepl('''#!/bin/bash
#$ -V                           # Inherit the submission environment
#$ -cwd                         # Start job in submission directory
#$ -N {{simulname}}             # Job Name
#$ -j y                         # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID        # Name of the output file (eg. myMPI.oJobID)
#$ -pe 12way 12                 # Requests 12 tasks/node, 12 cores total
#$ -q normal                    # Queue name normal
#$ -l h_rt=24:00:00             # Run time (hh:mm:ss) - 24.0 hours

# make sure to find the correct gromacs
export GMXLIB={{gmxtopdir}}

if [ f $prj.cpt ] ; then
  OPTCPT="-cpi $prj.cpt"
else
  OPTCPT=" "
fi

{{mdrun}} -maxh 23.9 -deffnm $prj $OPTCPT &

wait
''', d)



def main():
  global temp, pres, nequil, nsteps, nthreads
  global buildroot, gmxtopdir

  # make a working directory
  if not pres:
    prjname = "nvt"
    prjdir = "T%gdbl" % temp
    egrps = "11 0"
    fnxvg = "e.xvg"
  else:
    prjname = "npt"
    prjdir = "T%gP%gdbl" % (temp, pres)
    egrps = "11 21 0"
    fnxvg = "ev.xvg"
  if not os.path.isdir(prjdir):
    os.mkdir(prjdir)

  initdir = "initdbl"
  proggmx = "gmx_d"

  # copy files
  os.system("rm -f %s/*" % prjdir)
  shutil.copy(initdir + "/init.gro", prjdir + "/init.gro")
  shutil.copy(initdir + "/topol.top", prjdir + "/topol.top")

  # copy the NVT mdp
  s = mkmdp(initdir + "/nvt.mdp", temp)
  open(prjdir + "/nvt.mdp", "w").writelines(s)

  if pres:
    # copy the NPT mdp
    s = mkmdp(initdir + "/npt.mdp", temp, pres)
    open(prjdir + "/npt.mdp", "w").writelines(s)

  os.chdir(prjdir)

  envars = { "GMXLIB" : gmxtopdir }

  # NVT grompp
  changensteps("nvt.mdp", nequil)
  cmd = "%s/bin/%s grompp -f nvt.mdp -o nvt.tpr -c init.gro" % (
    buildroot, proggmx)
  zcom.shrun(cmd, envars = envars)

  # NVT mdrun
  cmd = "%s/bin/%s mdrun -v -deffnm nvt" % (
      buildroot, proggmx)
  if nthreads > 0:
    cmd += " -nt %d" % nthreads
  zcom.shrun(cmd, envars = envars)

  shutil.move("nvt.gro", "init.gro")

  if pres:
    # NPT grompp
    changensteps("npt.mdp", nequil)
    cmd = "%s/bin/%s grompp -f npt.mdp -o npt.tpr -c init.gro" % (
      buildroot, proggmx)
    zcom.shrun(cmd, envars = envars)

    # NPT mdrun
    cmd = "%s/bin/%s mdrun -v -deffnm npt" % (
        buildroot, proggmx)
    if nthreads > 0:
      cmd += " -nt %d" % nthreads
    zcom.shrun(cmd, envars = envars)

    shutil.move("npt.gro", "init.gro")

  # compute energy
  cmd = "echo %s | %s/bin/%s energy -f %s.edr -o %s" % (
      egrps, buildroot, proggmx, prjname, fnxvg)
  zcom.shrun(cmd, envars = envars)

  # remove unneeded files
  os.system("rm -f *.log *.cpt *.tpr *.xtc *.trr *mdout* \#*")

  # change the number of steps for production
  changensteps("%s.mdp" % prjname, nsteps)
  cmd = "%s/bin/%s grompp -f %s.mdp -o %s.tpr -c init.gro" % (
    buildroot, proggmx, prjname, prjname)
  zcom.shrun(cmd, envars = envars)

  os.chdir("..")

  #mklonestarpbs()



if __name__ == "__main__":
  doargs()
  main()
