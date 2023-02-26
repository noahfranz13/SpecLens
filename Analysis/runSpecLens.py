'''
Wrapper script for the SpecLens class to run
it from the command line.

By: Noah Franz
'''
import sys
import argparse
from SpecLens import SpecLens

def main():

    p = argparse.ArgumentParser()
    p.add_argument('-i', '--infile', help='input file path that includes columns for targetid, survey, program, and healpis', required=True)
    p.add_argument('--outdir', help='output directory, if none provided the default is the directory of the infile', default=None)
    p.add_argument('--specprod', help='spectroscopic production name, default is iron', default='iron')
    p.add_argument('--mp', help='how many cores to run with', default=1)
    p.add_argument('--overwrite', action='store_true', help='Overwrite any existing output files.')
    p.add_argument('--makeqa', action='store_true', help='Build QA in parallel.')
    args = p.parse_args()
    
    # initialize SpecLens object
    s = SpecLens(args.infile,
                 outdir=args.outdir,
                 specprod=args.specprod,
                 overwrite=args.overwrite,
                 mp=args.mp)

    s.separateLens() # separate all the lenses

    if args.makeqa:
        s.generateQA()

if __name__ == "__main__":
    sys.exit(main())
