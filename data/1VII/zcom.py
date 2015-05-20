#!/usr/bin/env python

''' commonly used python modules

    List of functions:
    * die_if(x)     die if `x' is true
    * runcmd(cmd)   run command and optionally capture the output
    * shrun(cmd)    run `cmd' with environmental variables
    * templrepl()   template replacement
    * mkpath()      make a directory
    * pathglob()    search files matching given patterns under a directory
    * argsglob()    glob command line arguments
'''

import os, sys, stat, shutil, getopt, re, subprocess, glob, random



def die_if(cond, message = ""):
  ''' raise an Exception, if `cond' is true '''
  if cond:
    print "fatal error:", message
    raise Exception



def runcmd(cmd, capture = 0, system = 0, verbose = 1):
  ''' run a system command and optionally capture standard output/error
  return a tuple (code, stdout, stderr)
  1. The latter two are None unless `capture' is set, and `system' is 0
  2. To use os.system() instead of subprocess, set `system' to 1
  3. If `verbose' is set, the command is echoed before executed '''

  pipe = None
  if capture: pipe = subprocess.PIPE

  # detect if `cmd' is a string or not
  if type(cmd) == str:
    cmd = cmd.split()

  if verbose >= 1:
    print "CMD:", ' '.join(cmd)
    if verbose >= 2: raw_input("proceed?")

  if system:
    ret = os.system( ' '.join(cmd) )
    oe = ["", ""]  # no stdout or stderr
  else:
    p = subprocess.Popen(cmd, stdout=pipe, stderr=pipe)
    oe = p.communicate()
    ret = p.returncode
  return (ret, oe[0], oe[1])



def shrun(cmd, capture = 0, verbose = 1,
          prefix = "#!/bin/bash", envars = {}, fnscript = None):
  ''' write a shell script to run the cmd `cmd'
      with environment variables in `envars' defined '''

  fn = "TMP" + str(random.randint(0, 99999)) +  ".sh"
  if fnscript: fn = fnscript
  s = prefix + "\n"
  for k in envars:
    s += "export " + k + "=" + envars[k] + "\n"
  s += cmd + "\n"
  open(fn, "w").write(s)
  os.chmod(fn, stat.S_IEXEC + stat.S_IREAD + stat.S_IWRITE)
  if verbose:
    print "CMD:", cmd
  ret = runcmd(os.path.abspath(fn), capture, verbose = 0)
  die_if (ret[0] != 0, "error occurred when running %s\nscript %s" % (cmd, fn))
  if not fnscript: os.remove(fn)
  return ret



def safebackup(fn0, ext = "", verbose = 1):
  ''' backup file `fn0' to a nonexisting name,
  lead by the original file name `fn.ext' + a numeric extension
  if `verbose' is >= 2, confirmation is needed after the backup '''

  fn = fn0 + ext
  i = 1
  if ext == "": fn += str(i)
  while os.path.exists(fn):
    fn = fn0 + ext + str(i)
    i += 1
  shutil.copy2(fn0, fn)
  if verbose:
    print fn0, "is backed up to", fn
    if verbose >= 2:
      raw_input("press Enter to proceed ...")



def templrepl(s, d, bra = "{{", ket = "}}", rex = False):
  ''' template replacement based on the dictionary `d'
      keywords are wrapped by `bra' and `ket'
      if `rex' is True, treat keys and values in the dictionary
      as regular expressions '''

  for k in d:
    key = k
    if not (key.startswith(bra) and key.endswith(ket)):
      key = bra + key + ket
    val = str( d[k] )
    if rex: # regular expression
      s = re.sub(key, val, s)
    else: # normal string
      s = s.replace(key, val)
  return s



def mkpath(path, make = 0, force = False):
  ''' return a new path as prefix + path, expand $HOME
  create the path as a directory, if `make' is true '''

  if type(path) == list:
    path = os.path.join(*path)
  path = os.path.realpath( os.path.expanduser(path) )

  if os.path.exists(path): # remove existing
    if not force: return path # skip an existing path
    shutil.rmtree(path)
  os.mkdir(path)
  print "making directory: ", path
  return path



def pglob(pat, files = True, dirs = False, links = False, dots = False):
  ''' list files that match `pat' '''

  ls = [ a for a in glob.glob(pat) ]
  if not files: # exclude files
    ls = [ a for a in ls if not os.path.isfile(a) ]
  if not dirs: # exclude directories
    ls = [ a for a in ls if not os.path.isdir(a) ]
  if not links: # exclude symbolic links
    ls = [ a for a in ls if not os.path.islink(a) ]
  if not dots: # exclude hidden files
    ls = [ a for a in ls if not a.startswith(".") ]
  # don't expand to real paths, because there may be aliases
  ## expand to real directires
  #ls = [ os.path.realpath(a) for a in ls ]
  # remove duplicate
  return sorted( list( set( ls ) ) )



def pathglob(pats, root = None, recur = False, files = True,
             dirs = False, links = False, dots = False):
  ''' find all files that match a list of patterns `pats' under `root'
      hidden directories `.git' are skipped unless `dots' is True '''

  if type(pats) == str: pats = [ pats ]

  if not root: root = os.getcwd()
  root = os.path.realpath( os.path.expanduser(root) )
  if not os.path.exists(root): return []

  ls = []
  if not recur: # list files under the current dir
    for pat in pats:
      ls += pglob(os.path.join(root, pat), files, dirs, links, dots)
  else: # recursively list all files under all subdirectories
    for r, subdirs, fns in os.walk(root):
      for pat in pats:
        ls += pglob(os.path.join(r, pat), files, dirs, links, dots)
      if not dots:
        # remove subdirectories that start with ".", such that
        # os.walk won't go into them further
        i = 0
        while i < len(subdirs):
          d = subdirs[i]
          if d.startswith("."):
            subdirs.remove(d)
          else: i += 1
  # remove duplicate
  return sorted( list( set(ls) ) )



def argsglob(args, defpats = "*", recur = False, links = False):
  ''' glob patterns in the command line arguments
      if no arguments are given, use the default `defpats',
      which is a string like "*.c *.txt" '''

  pats = defpats.split() # default pattern
  # if there are arguments, use the patterns in the argumens instead
  if len(args) > 0: # parse the pattern in each argument
    pats = [ a for pat in args for a in pat.split() ]
  return pathglob(pats, recur = recur, links = links)


