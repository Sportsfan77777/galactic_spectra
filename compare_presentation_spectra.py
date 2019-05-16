"""
Make spectra for presentations
"""

import math
import numpy as np
from scipy.interpolate import interp2d

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

def planck(wavelength, temperature):
    # Planck's Law
    wavelength_m = wavelength * 10**(-9) # convert nm to m
    return np.power(wavelength_m, -5.0) / np.expm1(h * c / k_B / wavelength_m / temperature)

###############################################################################

### Make plots ###

linewidth = 2
fontsize = 16
dpi = 250
tex = False
#grid_alpha = 0.4

def make_plot(show = False, num_blue = 1, num_yellow = 1, num_red = 1):
    fig = plot.figure(figsize = (7, 4))

    # Data
    wavelengths = np.linspace(100, 1000, 1000)

    spectrum_blue = num_blue * planck(wavelengths, 10000)
    spectrum_yellow = num_yellow * planck(wavelengths, 5700)
    spectrum_red = num_red * planck(wavelengths, 3500)
    spectrum = spectrum_blue + spectrum_yellow + spectrum_red

    # Plot spectra
    normalization = max(spectrum_red)

    x = wavelengths
    y = spectrum / normalization

    y_red = spectrum_red / normalization
    y_yellow = spectrum_yellow / normalization
    y_blue = spectrum_blue / normalization

    plot.plot(x, y, c = 'k', linewidth = linewidth + 1, zorder = 100)
    
    plot.plot(x, y_red, c = 'r', linewidth = linewidth)
    plot.plot(x, y_yellow, c = 'y', linewidth = linewidth)
    plot.plot(x, y_blue, c = 'b', linewidth = linewidth)

    # Grid Lines
    plot.grid(linestyle = "--")

    # Decorate with spectral colors
    alpha_rainbow = 0.03; num_colors = 500

    coordinates = np.linspace(400, 700, num_colors); y_region = np.array([10**(-6), 10000])
    visible_spectrum = np.zeros((num_colors, 2))
    visible_spectrum[:, 0] = coordinates; visible_spectrum[:, 1] = coordinates
    plot.pcolormesh(coordinates, y_region, np.transpose(visible_spectrum), cmap = 'nipy_spectral', alpha = 0.1)

    # Axes
    plot.xlim(100, 1000)
    plot.ylim(0, max(y))

    # Labels
    plot.xlabel(r"$\mathrm{Wavelength\ (nm)}$", fontsize = fontsize, usetex = tex)
    plot.ylabel(r"$\mathrm{Normalized\ Luminosity}$", fontsize = fontsize, usetex = tex)

    # Save + Show
    plot.savefig("presentation_spectra/comparison_test.png", dpi = dpi, bbox_inches = 'tight')
    if show:
        plot.show()

make_plot(show = True)

