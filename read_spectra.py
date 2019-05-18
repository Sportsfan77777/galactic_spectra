"""
Read real galaxy spectra data
"""

import os, sys
import glob

import math
import numpy as np
from scipy.interpolate import interp2d

from matplotlib import pyplot as plot

#fn = sys.argv[1]
#data = np.loadtxt(fn)
#wavelengths = data[:, 0] # lambda
#spectra = data[:, 1] # F_lambda

linewidth = 2
fontsize = 16
dpi = 200

### Helper Functions ###

### Plotting ###

def make_swire_plot():
    fig = plot.figure(figsize = (7, 4))

    fns = glob.glob("swire_library/*.sed")
    for fn in fns:
        data = np.loadtxt(fn)
        wavelengths = data[:, 0] / 10.0 # lambda (in nm)
        spectra = data[:, 1] # F_lambda

        name = fn.split("/")[-1].split(".")[0].split("_")[0]

        plot.plot(wavelengths, spectra, c = 'k', linewidth = linewidth)

        #plot.legend(loc = "upper right")

        #plot.xlim(10**(2), 10**(6))
        #plot.ylim(10**(-5), 10)

        plot.xlim(100, 1000)

        plot.xlabel("Wavelength (nm)", fontsize = fontsize)
        plot.ylabel("Normalized Luminosity", fontsize = fontsize)
        plot.title(name, fontsize = fontsize + 2)

        #plot.xscale('log')
        #plot.yscale('log')

        plot.savefig("spectra/swire-%s" % name)
        #plot.show()
        plot.cla()

def make_kc_plot():
    fig = plot.figure(figsize = (7, 4))

    fns = glob.glob("kc96/*.ascii")
    for fn in fns:
        data = np.loadtxt(fn)
        wavelengths = data[:, 0] / 10.0 # lambda (in nm)
        spectra = data[:, 1] # F_lambda

        normalization_index = np.searchsorted(wavelengths, 550)
        spectra /= spectra[normalization_index]

        name = fn.split("/")[-1].split(".")[0].split("_")[0]

        plot.plot(wavelengths, spectra, c = 'k', linewidth = linewidth)

        #plot.legend(loc = "upper right")

        #plot.xlim(10**(2), 10**(6))
        #plot.ylim(10**(-5), 10)

        #plot.xlim(100, 1000)

        plot.xlabel("Wavelength (nm)", fontsize = fontsize)
        plot.ylabel("Normalized Luminosity", fontsize = fontsize)
        plot.title(name, fontsize = fontsize + 2)

        #plot.xscale('log')
        #plot.yscale('log')

        plot.savefig("spectra/kc96-%s" % name, dpi = dpi)
        #plot.show()
        plot.cla()

#make_swire_plot()

make_kc_plot()
