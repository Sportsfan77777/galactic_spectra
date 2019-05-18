"""
Make spectra for presentations
"""

import os, sys

import math
import numpy as np
from scipy.interpolate import interp1d, interp2d

from matplotlib import pyplot as plot
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.colors

###############################################################################

# Constants
h = 6.626 * 10**(-34)
c = 3.0 * 10**(8)
k_B = 1.38 * 10**(-23)

###############################################################################

### Spectrum functions ###



###############################################################################

### Make plots ###

linewidth = 1
fontsize = 16
dpi = 250
tex = False
#grid_alpha = 0.4

def make_plot(show = False, log = False):
    fig = plot.figure(figsize = (7, 4))

    # Data
    fn = sys.argv[1]
    directory = fn.split("/")[0]
    name = fn.split("/")[-1].split(".")[0].split("_")[0]

    if directory == "kc96":
        # Kinsey + Calzetti 1996
        data = np.loadtxt(fn)
        wavelengths = data[:-15, 0] / 10.0 # lambda (in nm)
        spectrum = data[:-15, 1] # F_lambda

        normalization_index = np.searchsorted(wavelengths, 550)
        spectrum /= spectrum[normalization_index]
    else:
        # Swire Library
        data = np.loadtxt(fn)
        wavelengths = data[:, 0] / 10.0 # lambda (in nm)
        spectrum = data[:, 1] # F_lambda

    # Plot spectra
    x = wavelengths
    y = spectrum

    plot.plot(x, y, c = 'k', linewidth = linewidth, zorder = 100)

    # Grid Lines
    #plot.grid(linestyle = "--")

    # Decorate with spectral colors
    alpha_rainbow = 0.03; num_colors = 500

    coordinates = np.linspace(400, 700, num_colors); y_region = np.array([10**(-6), 100000])
    visible_spectrum = np.zeros((num_colors, 2))
    visible_spectrum[:, 0] = coordinates; visible_spectrum[:, 1] = coordinates
    plot.pcolormesh(coordinates, y_region, np.transpose(visible_spectrum), cmap = 'nipy_spectral', alpha = 0.1)

    # Axes
    plot.xlim(100, 1000)
    plot.ylim(10**(-6), 1.05 * max(spectrum))

    if log:
        plot.xlim(100, 10000)

    #plot.yscale('log')

    # Labels
    plot.xlabel(r"$\mathrm{Wavelength\ (nm)}$", fontsize = fontsize, usetex = tex)
    plot.ylabel(r"$\mathrm{Normalized\ Luminosity}$", fontsize = fontsize, usetex = tex)

    # Save + Show
    if log:
        plot.savefig("template_spectra/%s-%s-log.png" % (directory, name), dpi = dpi, bbox_inches = 'tight')
    else:
        plot.savefig("template_spectra/%s-%s.png" % (directory, name), dpi = dpi, bbox_inches = 'tight')

    if show:
        plot.show()

make_plot(show = True, log = True)

