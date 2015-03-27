#! /usr/bin/env python

"""
A program to provide a command line interface to the AegeanTools.fits_interp module

@author: Paul Hancock

Created:
3rd Feb 2015
"""

import logging
import argparse
import os
import sys
import numpy as np
from astropy.io import fits
from AegeanTools.fits_interp import compress, expand

__version__ = '1.9-0-gdae3d25'
__date__ = '2015-03-26'

# command line version of this program runs from here.
if __name__ == "__main__":

    epilog = ''
    parser = argparse.ArgumentParser(epilog=epilog, prefix_chars='-')

    # tools for shrinking files
    group1 = parser.add_argument_group("Shrinking and expanding files")
    group1.add_argument('infile', type=str,
                        help="input filename")
    group1.add_argument('-o', dest='outfile', action='store',
                        default=None, type=str, metavar='OutputFile',
                        help='output filename')
    group1.add_argument('-f', dest='factor', action='store',
                        default=None, type=int, metavar='factor',
                        help='reduction factor')
    group1.add_argument('-x', dest='expand', action='store_true',
                        default=False,
                        help='Operation is expand instead of compress.')
    group1.add_argument('-i', dest='mode', choices=['linear','nearest','cubic'],
                        default='linear',
                        help='Interpolation method')
    group1.add_argument('-m', dest='maskfile',action='store',
                        default=None, type=str, metavar='MaskFile',
                        help="File to use for masking pixels.")
    # TODO: move these to be in a different group. (The same as help).
    group1.add_argument('--debug', dest='debug', action='store_true',
                        default=False,
                        help='Debug output')
    group1.add_argument('--version', action='version', version='%(prog)s '+__version__+"-({0})".format(__date__))

    results = parser.parse_args()

    logging_level = logging.DEBUG if results.debug else logging.INFO
    logging.basicConfig(level=logging_level, format="%(process)d:%(levelname)s %(message)s")
    logging.info("This is SR6 {0}-({1})".format(__version__,__date__))


    if not os.path.exists(results.infile):
        logging.error("{0} does not exist".format(results.infile))
        sys.exit()

    if results.expand:
        if results.maskfile and os.path.exists(results.maskfile):
            maskdata = fits.open(results.maskfile)[0].data
            mask = np.where(np.isnan(maskdata))
            hdulist = expand(results.infile, method=results.mode)
            hdulist[0].data[mask]=np.nan
            hdulist.writeto(results.outfile,clobber=True)
            logging.info("Wrote masked file: {0}".format(results.outfile))
        elif results.maskfile is None:
            expand(results.infile, results.outfile, method=results.mode)
        else:
            logging.error("Can't find {0}".format(results.maskfile))

    else:
        compress(results.infile, results.factor, results.outfile)
