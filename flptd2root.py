#!/usr/bin/env python
import os
import argparse
import subprocess

def getConversionScript():
    """ derive path to conversion script from this
    """
    thisScript = os.path.realpath(__file__)
    return os.path.join(os.path.dirname(thisScript), "p2r.conf")


def flbrio2root():
    """implementation of Brio->Root conversion
    """
    parser = argparse.ArgumentParser(description="Convert flreconstruct BRIO format to ROOT")
    parser.add_argument("-i", type=str, metavar="<infile>", required=True, help="input BRIO file")
    parser.add_argument("-o", type=str, metavar="<outfile>", required=True, help="output ROOT file")
    args = parser.parse_args()

    try:
        subprocess.call(["flreconstruct", "-i", args.i, "-p", getConversionScript()])
        os.rename("testptd.root", args.o)
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            print("error: could not locate flreconstruct program")
        else:
            # Something else went wrong while trying to run `flreconstruct`
            raise


if __name__ == '__main__':
    flbrio2root()
