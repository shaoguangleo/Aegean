"""
Little script to compare two fits image files and determine whether they are similar to a certain precision.
"""
#! /usr/bin/env python

import sys
import numpy as np
from astropy.io import fits

DEFAULT_PRECISION = 0.00000001

def get_data(filename):
	"""
	Returns an ndarray containing the data from the file specified.
	"""
	# Figure out how many axes are in the datafile, and the shape.
	header = fits.getheader(filename)
	x_max = header['NAXIS2']
	y_max = header['NAXIS1']
	NAXIS = header["NAXIS"]

	# It seems that I cannot memmap the same file multiple times without errors
	with fits.open(filename, memmap=False) as a:
		if NAXIS == 2:
			data = a[0].section[0:x_max, 0:y_max]
		elif NAXIS == 3:
			data = a[0].section[0, 0:x_max, 0:y_max]
		elif NAXIS == 4:
			data = a[0].section[0, 0, 0:x_max, 0:y_max]
		else:
			print "Too many NAXIS for me {0}".format(NAXIS)
			print "Fix file {0} to be more sane".format(filename)
			sys.exit(1)
	return data

if __name__ == "__main__":
	if len(sys.argv) < 0:
		print "Usage: compare_fits.py <first.fits> <second.fits>"
	first = get_data(sys.argv[1])
	second = get_data(sys.argv[2])
	if first.shape != second.shape:
		print "Files differ: shapes are different.", first.shape, second.shape
		sys.exit(1)
	# First make a mask of the NaN entries, since we can't meaningfully subtract these.
	mask = np.isfinite(first)
	if (np.isfinite(second) != mask).any():
		print "Files differ: different arrangement of non-finite values (NaN, Inf etc.)"
		sys.exit(1)
	# Now we can assume that both files have the same mask, so use it to do our subtraction.
	masked_first = first[mask]
	matrix_diff = masked_first[abs(masked_first - second[mask]) > 0.00000001]
	# matrix_diff = first[abs(first - second) > DEFAULT_PRECISION]
	if matrix_diff.size:
		total_entries = first.shape[0]*first.shape[1]
		print "Files differ: a total of [%s/%s] (%5.1f%%) entries are outside the allowed precision." % (matrix_diff.size,
			total_entries, ((matrix_diff.size * 100.0) / total_entries))
		sys.exit(1)
	# Files are close enough. :)
	sys.exit(0)