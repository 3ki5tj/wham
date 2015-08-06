#!/usr/bin/env python


import os, sys, codecs



rmcomments = False
rmspaces = False



def mkdic(fntrans):
  fp = codecs.open(fntrans, "r", "utf-8")
  stab = fp.read()
  if stab[0] == u'\ufeff':
    stab = stab[1:]
  fp.close()
  dic = {}
  arr = stab.split('\n')
  for i in range(len(arr)):
    ln = arr[i].strip()
    if ln.find('\t') < 0:
      continue
    a = ln.split('\t')
    dic[a[0]] = a[1]
  return dic



def translate(src, dic):
  ''' make a translation '''
  if dic == None:
    return

  for i in range(len(src)):
    if src[i].strip() == "":
      continue
    for k in dic:
      src[i] = src[i].replace(k, dic[k])



def loadsrc(fn):
  ''' load a script or style sheet '''

  src = open(fn).readlines()
  for i in range(len(src)):
    ln = src[i].strip()
    if ln == "":
      src[i] = "" # remove blank lines
    elif rmcomments:
      if ln.startswith("/*"):
        # remove block comments
        for j in range(i, len(src)):
          if src[j].strip().endswith("*/"):
            break
        for k in range(i, j+1):
          src[k] = ""
      else:
        # remove line comments
        p = src[i].rfind("//")
        if p >= 0:
          src[i] = src[i][:p] + "\n"

  # remove spaces
  if rmspaces:
    for i in range(len(src)):
      ln = src[i].strip()
      if src[i].startswith("}") or ln == '"use strict";':
        ln += "\n"
      src[i] = ln

  return ''.join(src)


def htmlpack(fn, fntrans = None):
  ''' embed scripts and stylesheets into the HTML file `fn` '''

  s = open(fn).readlines()
  n = len(s)
  for i in range(n):
    ln = s[i].strip()
    if ln.startswith("<script"):
      p = ln.find('src="')
      p = ln.find('"', p)
      q = ln.find('"', p+1)
      fnsrc = ln[p+1:q]
      s[i] = ('<!-- %s -->\n' % fnsrc) + '<script>\n' + loadsrc(fnsrc) + '</script>\n'
    elif ln.startswith('<link rel='):
      p = ln.find('href="')
      p = ln.find('"', p)
      q = ln.find('"', p+1)
      fnsrc = ln[p+1:q]
      s[i] = ('<!-- %s -->\n' % fnsrc) + '<style>\n' + loadsrc(fnsrc) + '</style>\n'

  arr = os.path.splitext(fn)
  fnout = arr[0] + "_pack" + arr[1]
  open(fnout, "w").writelines(s)
  print "saved to " + fnout

  if fntrans and os.path.exists(fntrans):
    # translation
    dic = mkdic(fntrans)
    translate(s, dic)
    s = ''.join(s)
    fnout = arr[0] + "_pack_cn" + arr[1]
    fp = codecs.open(fnout, 'w', 'utf-8')
    fp.write(u'\ufeff')
    fp.write(s)
    fp.close()
    print "saved to " + fnout



if __name__ == "__main__":
  fn = "ljmix.html"
  if len(sys.argv) > 1:
    fn = sys.argv[1]
  print fn
  htmlpack(fn)
