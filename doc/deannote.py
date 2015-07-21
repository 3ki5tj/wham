#!/usr/bin/env python



r''' remove annotations in the tex file

The applicable macros are the following

\usepackage[usenames,dvipsnames]{xcolor}

...

\newcommand{\repl}[2]{{\color{gray} [#1] }{\color{blue} #2}}
\newcommand{\add}[1]{{\color{blue} #1}}
\newcommand{\del}[1]{{\color{gray} [#1]}}
\newcommand{\note}[1]{{\color{OliveGreen}\small [\textbf{Comment.} #1]}}

'''



import sys, os, getopt, shutil, re



fninp = None
fnout = None
rmcom = False
verbose = 0



def usage():
  ''' print usage and die '''

  print """  Usage:
    %s [OPTIONS]""" % sys.argv[0]

  print """
  remove annotations in the tex file

  OPTIONS:
    -i              set the input file
    -o              set the output file
    -c              remove comments
    -v              be more verbose
    --verbose=      set verbosity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "i:o:cv",
        [ "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fninp, fnout, rmcom, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-i",):
      fninp = a
    elif o in ("-o",):
      fnout = a
    elif o in ("-c",):
      rmcom = True
    elif o in ("-h", "--help"):
      usage()

  if len(args):
    fninp = args[0]

  if len(args) >= 2:
    fnout = args[1]



def findbraceend(s, p):
  ''' find the ending brace '''
  level = 1
  for i in range(p, len(s)):
    if s[i] == "{":
      level += 1
    elif s[i] == "}":
      level -= 1
      if level <= 0:
        break
  return i



def eatprev(s, p):
  ''' eat the previous newline '''
  if (p >= 1 and s[p - 1] == '\n'
      and p < len(s) - 1 and s[p + 1] == '\n'):
    p -= 1
  return p



def eatnext(s, p):
  ''' eat the next newline '''
  if (p >= 1 and s[p - 1] == '\n'
      and p < len(s) - 1 and s[p + 1] == '\n'):
    p += 1
  return p + 1



def add(tag, s):
  while 1:
    p = s.find(tag)
    if p < 0:
      break
    q = findbraceend(s, p + len(tag))
    if q < 0:
      print "cannot find brace ending for the added text! %s, %s" % (p, tag)
      break
    if verbose:
      print "Adding [%s]" % repr(s[p+len(tag):q].strip())
      if verbose >= 2: raw_input()
    q1 = q + 1
    s = s[:p] + s[p+len(tag):q] + s[q1:]
  return s



def delete(tag, s):
  while 1:
    p = s.find(tag)
    if p < 0:
      break
    q = findbraceend(s, p + len(tag))
    if q < 0:
      print "cannot find brace ending for the deleted text! %s, %s" % (p, tag)
      break
    p0, q0 = p, q
    q1 = eatnext(s, q)
    if verbose:
      print "Removing [%s]" % repr(s[p0+len(tag):q].strip())
      if verbose >= 2: raw_input()
    s = s[:p] + s[q1:]
  return s



def repl(tag, s):
  while 1:
    p = s.find(tag)
    if p < 0:
      break
    q = findbraceend(s, p + len(tag))
    if q < 0:
      print "cannot find brace ending for the replaced text! %s, %s" % (p, tag)
      break

    # find the new text
    pp = s.find("{", q+1)
    qq = findbraceend(s, pp+1)
    if qq < 0:
      print "cannot find brace ending for the candidate text! %s" % (pp)
      break

    p0 = p
    pp1 = eatnext(s, pp)
    qq1 = eatnext(s, qq)
    if verbose:
      print "Replacing [%s] by [%s]" % (repr(s[p0+len(tag):q].strip()), repr(s[pp1:qq].strip()))
      if verbose >= 2: raw_input()
    s = s[:p] + s[pp1:qq] + s[qq1:]
  return s



def rmcomments(s):
  arr = s.split('\n')
  arr2 = []
  for i in range(len(arr)):
    if not arr[i].strip().startswith("%"):
      arr2 += [ arr[i], ]
  return '\n'.join(arr2)



def main():
  global fninp, fnout

  s = open(fninp).read()

  s = delete(r"\note{", s)
  s = delete(r"\del{", s)
  s = add(r"\add{", s)
  s = repl(r"\repl{", s)

  # remove trailing spaces of each line
  arr = s.split('\n')
  for i in range(len(arr)):
    arr[i] = arr[i].rstrip()
  s = '\n'.join(arr)

  if rmcom:
    s = rmcomments(s)

  if not fnout:
    arr = os.path.splitext(fninp)
    fnout = arr[0] + "_clean" + arr[1]

  print "saving results to %s" % fnout
  open(fnout, "w").write(s)



if __name__ == "__main__":
  doargs()
  main()

