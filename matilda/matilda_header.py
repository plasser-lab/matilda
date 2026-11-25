from __future__ import print_function, division
import os

width=80

def print_header(*args, **kwargs):
    print((ret_header(*args, **kwargs)))

def ret_header(title=None, ioptions=None, cfile=None, ver='1.0'):
    hstr  = width*'=' + '\n'

    hstr += addlinec("MATILDA %s"%ver)
    hstr += addlinec()
    hstr += addlinec("Authors: F. Plasser, S. Ghosh, U. Rehman, G. Woelfle Conway")

    if not title==None:
        hstr += width*'-' + '\n'
        hstr += addlinec(title)

    hstr += width*'=' + '\n'

    return hstr

def addlinec(line=""):
    return "|" + line.center(width-2) + "|\n"

def addlinel(line="", lpad=5):
    return "|" + lpad*' ' + line.ljust(width-2-lpad) + "|\n"
