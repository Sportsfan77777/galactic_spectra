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

### Plotting functions ###

def wavelength_to_rgb(wavelength, gamma = 1):
    # https://stackoverflow.com/questions/44959955/matplotlib-color-under-curve-based-on-spectral-color
    ''' Taken from http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
    This converts a given wavelength of light to an approximate RGB color value. 
    The wavelength must be given in nanometers in the range from 380 nm through 750 nm (789 THz through 400 THz).

    Based on code by Dan Bruton (http://www.physics.sfasu.edu/astro/color/spectra.html)
    Additionally alpha value set to 0.5 outside range'''

    if wavelength >= 380 and wavelength <= 750:
        A = 1.0
    else:
        A = 1.0
    if wavelength < 380:
        wavelength = 380.
    if wavelength >750:
        wavelength = 750.
    if wavelength >= 380 and wavelength <= 440:
        attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
        R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ** gamma
        G = 0.0
        B = (1.0 * attenuation) ** gamma
    elif wavelength >= 440 and wavelength <= 490:
        R = 0.0
        G = ((wavelength - 440) / (490 - 440)) ** gamma
        B = 1.0
    elif wavelength >= 490 and wavelength <= 510:
        R = 0.0
        G = 1.0
        B = (-(wavelength - 510) / (510 - 490)) ** gamma
    elif wavelength >= 510 and wavelength <= 580:
        R = ((wavelength - 510) / (580 - 510)) ** gamma
        G = 1.0
        B = 0.0
    elif wavelength >= 580 and wavelength <= 645:
        R = 1.0
        G = (-(wavelength - 645) / (645 - 580)) ** gamma
        B = 0.0
    elif wavelength >= 645 and wavelength <= 750:
        attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
        R = (1.0 * attenuation) ** gamma
        G = 0.0
        B = 0.0
    else:
        R = 0.0
        G = 0.0
        B = 0.0
    return (R, G, B, A)

start_wavelength = 379; end_wavelength = 751
clim = (start_wavelength, end_wavelength)
norm = plot.Normalize(*clim)
spectrum_wavelengths = np.linspace(start_wavelength, end_wavelength, 1000)

color_list = list(zip(norm(spectrum_wavelengths),[wavelength_to_rgb(wavelength_i) for wavelength_i in spectrum_wavelengths]))
spectral_map = LinearSegmentedColormap.from_list("spectrum", color_list)

###############################################################################

### Make plots ###

linewidth = 2
fontsize = 17
dpi = 250
tex = False
dark_alpha = 1.0 # 0.75 hides most of the color

def make_plot(show = False):
    fig = plot.figure(figsize = (8, 4))

    # Data
    wavelengths = np.linspace(100, 1000, 1000)
    spectrum = planck(wavelengths, 15700)

    # Spectral grid
    luminosities = np.linspace(0, 2, 10)
    x_grid, y_grid = np.meshgrid(wavelengths, luminosities)
    extent = (wavelengths[0], wavelengths[-1], luminosities[0], luminosities[-1])


    # Plot spectra
    x = wavelengths
    y = spectrum / max(spectrum)

    uv_end = np.searchsorted(wavelengths, start_wavelength)
    ir_start = np.searchsorted(wavelengths, end_wavelength)

    plot.plot(x, y, c = 'k', linewidth = linewidth)

    # Grid Lines
    plot.grid(linestyle = '--')

    # Decorate with spectral colors
    plot.imshow(x_grid, extent = extent, clim = clim, cmap = spectral_map, aspect = 'auto')

    upper_curve = 2
    plot.fill_between(x, y, upper_curve, color = 'w')

    plot.fill_between(x[:uv_end], y[:uv_end], color = 'k', alpha = dark_alpha)
    plot.fill_between(x[ir_start:], y[ir_start:], color = 'k', alpha = dark_alpha)

    # Axes
    plot.xlim(100, 1000)
    plot.ylim(0, 1)

    # Labels
    plot.xlabel(r"$\mathrm{Wavelength\ (nm)}$", fontsize = fontsize, usetex = tex)
    plot.ylabel(r"$\mathrm{Normalized\ Luminosity}$", fontsize = fontsize, usetex = tex)

    # Save + Show
    plot.savefig("presentation_spectra/test.png", dpi = dpi, bbox_inches = 'tight')
    if show:
        plot.show()

make_plot(show = True)

