"""
Make spectra for presentations
"""

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

def planck(wavelength, temperature):
    # Planck's Law
    wavelength_m = wavelength * 10**(-9) # convert nm to m
    return np.power(wavelength_m, -5.0) / np.expm1(h * c / k_B / wavelength_m / temperature)

def luminosity(wavelength, temperature):
    # Luminosity: 4 * pi * R^2 * \sigma * T^4, where \sigma = 1 and R is in solar radii
    def interpolate_radius(temperature):
        # Using reference values from https://en.wikipedia.org/wiki/Main_sequence
        temperatures = np.array([2660, 3120, 3920, 4410, 5240, 5610, 5780, 5920, 6540, 7240, 8620, 10800, 16400, 30000, 38000])
        radii = np.array([0.13, 0.32, 0.63, 0.74, 0.85, 0.93, 1.0, 1.05, 1.2, 1.3, 1.7, 2.5, 3.8, 7.4, 18.0])

        interpolator = interp1d(temperatures, radii)
        radius = interpolator(temperature)

        return radius

    radius = interpolate_radius(temperature)
    return 4.0 * np.pi * np.power(radius, 2) * planck(wavelength, temperature)

###############################################################################

### Make plots ###

linewidth = 2
fontsize = 16
dpi = 250
tex = False
#grid_alpha = 0.4

def make_plot(show = False, num_blue = 1, num_yellow = 1, num_red = 1, ylim = "red"):
    fig = plot.figure(figsize = (7, 4))

    # Data
    wavelengths = np.linspace(100, 1000, 1000)

    T_blue = 10500
    T_yellow = 5700
    T_red = 3200

    one_blue = luminosity(wavelengths, T_blue)
    one_yellow = luminosity(wavelengths, T_yellow)
    one_red = luminosity(wavelengths, T_red)

    spectrum_blue = num_blue * one_blue
    spectrum_yellow = num_yellow * one_yellow
    spectrum_red = num_red * one_red
    spectrum = spectrum_blue + spectrum_yellow + spectrum_red

    # Plot spectra
    normalization = max(one_red)

    x = wavelengths
    y = spectrum / normalization

    y_red = spectrum_red / normalization
    y_yellow = spectrum_yellow / normalization
    y_blue = spectrum_blue / normalization

    plot.plot(x, y, c = 'k', linewidth = linewidth + 1, zorder = 100, label = "total")

    if num_red > 0:
        plot.plot(x, y_red, c = 'r', linewidth = linewidth)
    if num_yellow > 0:
        plot.plot(x, y_yellow, c = 'y', linewidth = linewidth)
    if num_blue > 0:
        plot.plot(x, y_blue, c = 'b', linewidth = linewidth)

    # Grid Lines
    plot.grid(linestyle = "--")

    # Decorate with spectral colors
    alpha_rainbow = 0.03; num_colors = 500

    coordinates = np.linspace(400, 700, num_colors); y_region = np.array([10**(-6), 100000])
    visible_spectrum = np.zeros((num_colors, 2))
    visible_spectrum[:, 0] = coordinates; visible_spectrum[:, 1] = coordinates
    plot.pcolormesh(coordinates, y_region, np.transpose(visible_spectrum), cmap = 'nipy_spectral', alpha = 0.1)

    # Axes
    plot.xlim(100, 1000)
    if ylim == "red":
        plot.ylim(10**(-6), max(y_red))
    elif ylim == "yellow":
        plot.ylim(10**(-6), max(y_yellow))
    else:
        plot.ylim(10**(-6), max(y))

    #plot.yscale('log')

    # Labels
    plot.xlabel(r"$\mathrm{Wavelength\ (nm)}$", fontsize = fontsize, usetex = tex)
    plot.ylabel(r"$\mathrm{Normalized\ Luminosity}$", fontsize = fontsize, usetex = tex)

    plot.legend()

    # Save + Show
    plot.savefig("presentation_spectra/comparison_test_%dB(%d)-%dY(%d)-%dR(%d)-limit%s.png" % (num_blue, T_blue, num_yellow, T_yellow, num_red, T_red, ylim), dpi = dpi, bbox_inches = 'tight')
    if show:
        plot.show()

num_red = 1
num_yellow = 0
num_blue = 0
ylim = "all"
make_plot(show = True, num_red = num_red, num_yellow = num_yellow, num_blue = num_blue, ylim = ylim)

