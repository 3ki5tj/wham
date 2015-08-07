#!/usr/bin/env python



r''' remove url and doi of the .bib file 

'''



import sys, os, getopt, shutil, re



fninp = None
fnout = None
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
    -v              be more verbose
    --verbose=      set verbosity
    -h, --help      help
  """
  exit(1)



def doargs():
  ''' handle input arguments '''
  try:
    opts, args = getopt.gnu_getopt(sys.argv[1:],
        "i:o:v",
        [ "help", "verbose=",
        ] )
  except getopt.GetoptError, err:
    print str(err)
    usage()

  global fninp, fnout, verbose

  for o, a in opts:
    if o in ("-v",):
      verbose += 1  # such that -vv gives verbose = 2
    elif o in ("--verbose",):
      verbose = int(a)
    elif o in ("-i",):
      fninp = a
    elif o in ("-o",):
      fnout = a
    elif o in ("-h", "--help"):
      usage()

  if len(args):
    fninp = args[0]

  if len(args) >= 2:
    fnout = args[1]



def delete(s, tag):
  ntag = len(tag)
  for i in range(len(s)):
    ln = s[i].strip()
    tag1 = ln[:ntag].lower()
    if tag1 == tag:
      if not ln.endswith(","):
        j = i - 1
        while j >= 0:
          if s[j].strip() != "":
            break
          j -= 1
        if j >= 0 and s[j].strip().endswith(","):
          s[j] = s[j].rstrip()[:-1] + "\n"
      s[i] = ""
          
  return s




def main():
  global fninp, fnout

  s = open(fninp).readlines()

  s = delete(s, "url")
  s = delete(s, "doi")

  if not fnout:
    arr = os.path.splitext(fninp)
    fnout = arr[0] + "_nourl" + arr[1]

  print "saving results to %s" % fnout
  open(fnout, "w").writelines(s)



if __name__ == "__main__":
  doargs()
  main()

