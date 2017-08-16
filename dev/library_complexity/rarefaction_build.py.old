#!/usr/bin/python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Description:	generate rarefaction curve.
#				From arXiv:1511.07428 and based on Dr. Silvano Garnerone's code.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
from joblib import Parallel, delayed
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# PARAMETERS ===================================================================

# Add script description
parser = argparse.ArgumentParser(
	description = 'Produces rarefaction curve.'
)

# Add params
parser.add_argument('dirpath', metavar = 'dirpath', type = str, nargs = 1,
	help = 'Experiment condition directory, contains UMI counts.')
parser.add_argument('experiment', metavar = 'expID', type = str, nargs = 1,
	help = 'Experiment ID (e.g., TK26).')
parser.add_argument('condition', metavar = 'cond', type = str, nargs = 1,
	help = 'Condition folder name (e.g., 400U2h).')

# Add flags
parser.add_argument('--step', metavar = 'step', type = float, nargs = 1,
	default = [0.1],
	help = 'Step size on the x-axis of the curve. Default: 0.1')
parser.add_argument('--niter', metavar = 'n', type = int, nargs = 1,
	default = [1000],
	help = 'Number of samples for the smoothing distro. Default: 1000')
parser.add_argument('--threads', metavar = 'threads', type = int, nargs = 1,
	default = [1], help = 'Number of threads for parallelization.')
parser.add_argument('-l',
	action = 'store_const', const = True, default = False,
	help = 'Positions were not converted to cutsite using pos2cs.')

# Parse arguments
args = parser.parse_args()

# Retrieve arguments
dirpath = args.dirpath[0]
experiment = args.experiment[0]
condition = args.condition[0]
step = args.step[0]
sd_samples = args.niter[0]
threads = args.threads[0]
atcs = not args.l

# Select file
if atcs:
	fname = "UMIpos.unique.atcs.txt"
else:
	fname = "UMIpos.unique.txt"

# FUNCTIONS ====================================================================

def smoothedGT(Phi, n, t, sd_samples):
	'''Compute the smoothed Good-Toulmin estimator.
	Based on arXiv:1511.07428 and Dr. Silvano Garnerone's code.
	Using binomial smoothing with q = 2/(2+t).


	Args:
		Phi (np.ndarray): list of prevalences.
		n (int): sample size.
		t (float): fold-changes of the new sample size wrt n.
		sd_samples (int): number of samples from the smoothing distribution.

	Returns:
		float: the smoothed Good-Toulmin estimator.
	'''

	# Scenario 1
	if t <= 1:
		coeff = np.array([-(-t)**(i + 1) for i in range(Phi.shape[0])])
		U_GT = sum(coeff * Phi)
		return np.mean(U_GT)

	# Scenario 2
	if t > 1:
		k = math.ceil( 0.5 * math.log(1. * n * t**(2 / (t - 1)), 3))
		q = 2. / (2 + t)
		L_list = np.random.binomial(k, q, size = (sd_samples, 1))
		U_list = [poisson_smoothing(t, L, Phi) for L in L_list]
		return np.mean(U_list)

def poisson_smoothing(t, L, Phi):
	'''Perform Poissonian smoothing.

	Args:
		t (float): fold-changes of the new sample size wrt n.
		L (int): binomial random variable (trimming).
		Phi (np.ndarray): list of prevalences.

	Returns:
		float: Good-Toulmin estimator.
	'''

	coeff = np.array([-(-t)**(i + 1) for i in range(L[0])])
	
	# Fix arrays shape
	if 0 < coeff.shape[0] - Phi.shape[0]:
		Phi = np.lib.pad(Phi, (0, coeff.shape[0] - Phi.shape[0]),
			'constant', constant_values = 0)

	return(sum(coeff * Phi[:L[0]]))

def build_curve(hist, n, x, sd_samples):
	'''Calculate the curve point.

	Args:
		Phi (np.ndarray): list of prevalences.
		n (int): sample size.
		t (float): fold-changes of the new sample size wrt n.
		sd_samples (int): number of samples from the smoothing distribution.

	Returns:
		float: the smoothed Good-Toulmin estimator.w
	'''
	if x == 0:
		return hist.sum()
	if x > 0:
		return hist.sum() + smoothedGT(hist, n, x, sd_samples)

# RUN ==========================================================================

# Read -------------------------------------------------------------------------

# Add trailing slash
if not "/" == dirpath[-1]:
	dirpath += "/"

# Read UMI pos
data = pd.read_csv(dirpath + fname, sep = "\t",
	names = ["chr", "starts", "seq", "counts"])

# Pre-processing ----------------------------------------------------------------

# Get occurrences of duplicate counts
counts = []
[counts.extend(n.split(" ")) for n in data['counts']]
counts = np.array(counts)
hist_partial = np.unique(counts, return_counts = True)

# Add missing numbers (# duplicates corresponds to the index)
hist = np.zeros(hist_partial[0].astype('int').max())
hist[hist_partial[0].astype('int') - 1] = hist_partial[1]

# Parameters
n = int(hist.sum())			# Number of duplicates
t = int(math.log(n, 10))	# Maximum fold-change for the curve (log(n))

# Prepare accumulation curve ---------------------------------------------------

# Curve coordinates
xs = np.arange(0.0, t, step)
ys = []

# Calculate ordinate coordinates
ys = Parallel(n_jobs = threads)(
			delayed(build_curve)(hist, n, x, sd_samples)
			for x in xs)
ys = np.array(ys)

# Output -----------------------------------------------------------------------

# Output dataset
pd.DataFrame(np.transpose(np.array([xs, ys])), columns = ['x', 'y']).to_csv(
	dirpath + experiment + '.' + condition + '.rarefaction.tsv',
	'\t', index = False)

# Output plot
plt.plot(xs, ys)
plt.suptitle('[' + experiment + '] Good-Toulmin based rarefaction curve.',
	fontsize=14, fontweight='bold')
plt.title('n.iter=%d ; step=%f' % (sd_samples, step))
plt.xlabel('Duplication level (as number of UMI occurrences).')
plt.ylabel('Estimated number of unseen amplification duplicates.')
pp = PdfPages(dirpath + experiment + '.' + condition + '.rarefaction.pdf')
plt.savefig(pp, format = 'pdf')
pp.close()

# END ==========================================================================

################################################################################
