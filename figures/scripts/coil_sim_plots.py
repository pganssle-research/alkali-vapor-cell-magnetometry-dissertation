"""
Generate pulse coil simulation plots

@author Paul J. Ganssle
"""
from __future__ import division

from scipy.io import loadmat
from numpy import *;
from os.path import join as pathjoin

from matplotlib import rcParams, rc
from matplotlib.pyplot import *;

from matplotlib.patches import Rectangle, FancyArrowPatch, Arc
import matplotlib.patheffects as PathEffects
from matplotlib.lines import Line2D


from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import theme.colors as thc
import theme.figure_settings as tfs

from matplotlib.colors import ColorConverter

import argparse

from general_util import maxf, minf, find_where

parser = argparse.ArgumentParser(description="Generate pulse sequence acquisitions.")
parser.add_argument('-s', '--save', help='Save the figures', action='store_true')
parser.add_argument('-c', '--close', help='Close the figures when done', action='store_true')
args = parser.parse_args()

SaveFigs = args.save
CloseFigs = args.close

tfs.setup_environment()

loc = '../coils'
make_path = lambda name: pathjoin(loc, name)
formats = ['png', 'eps']

figsize=(13, 3)
dpi = 150

# Load data
fname = 'data_sets/2012-08-16-GradCoilImage.mat'
o = loadmat(fname)
