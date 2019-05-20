"""
Make spectra for presentations

Note: Reference stellar spectra: http://vizier.u-strasbg.fr/viz-bin/VizieR-4
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

# Variables
T_sun = 5780

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

def solar_spectrum_interpolator():
    data = np.loadtxt("solar_spectrum.txt")
    wavelengths_sun = data[:, 0]
    spectrum_sun = data[:, 1]

    normalization = luminosity(1000, T_sun) / spectrum_sun[-1] # normalize to \lambda = 1000 nm, where emissivity = 1 roughly
    spectrum_sun *= normalization

    interpolator = interp1d(wavelengths_sun, spectrum_sun)
    return interpolator
solar_spectrum = solar_spectrum_interpolator()

###############################################################################

### Make plots ###

linewidth = 2
fontsize = 16
dpi = 250
tex = False
#grid_alpha = 0.4

def make_plot(show = False, comparison = True):
    fig = plot.figure(figsize = (7, 4))

    # Data
    wavelengths = np.linspace(120, 1000, 1000)

    solar_blackbody = luminosity(wavelengths, T_sun)
    real_solar_spectrum = solar_spectrum(wavelengths)

    # Plot spectra
    normalization = max(solar_blackbody)

    solar_blackbody /= normalization
    real_solar_spectrum /= normalization

    x = wavelengths
    y1 = solar_blackbody
    y2 = real_solar_spectrum

    plot.plot(x, y1, c = 'y', linewidth = linewidth, zorder = 100, label = "perfect blackbody")
    if comparison:
        plot.plot(x, y2, c = 'k', linewidth = linewidth - 1, zorder = 100, label = "real spectrum")

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
    plot.ylim(0, 1.2)

    #plot.yscale('log')

    # Labels
    plot.xlabel(r"$\mathrm{Wavelength\ (nm)}$", fontsize = fontsize, usetex = tex)
    plot.ylabel(r"$\mathrm{Normalized\ Luminosity}$", fontsize = fontsize, usetex = tex)

    plot.legend()

    # Save + Show
    if comparison:
        plot.savefig("presentation_spectra/real_vs_blackbody_%d.png" % T_sun, dpi = dpi, bbox_inches = 'tight')
    else:
        plot.savefig("presentation_spectra/blackbody.png", dpi = dpi, bbox_inches = 'tight')

    if show:
        plot.show()

make_plot(show = True, comparison = False)

