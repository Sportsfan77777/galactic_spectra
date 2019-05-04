"""
Plot dust emission (from Dale, Helou, et al. 2014???)
"""

import math
import numpy as np
from scipy.interpolate import interp2d
from matplotlib import pyplot as plot

num_files = 64

dust_emission_data = np.zeros((num_files + 1, 215))
dust_emission_data[:-1, 0] = -35
dust_emission_data[-1, :] = -50

wavelengths = np.ones(215)
wavelengths[0] = 0.05

def read_data():
    for i in range(64):
        fn = "dust_emission/spec_%02d.dat" % (i + 1)
        data_i = np.loadtxt(fn)
        wavelengths[1:] = data_i[:, 0]
        dust_emission_data[i, 1:] = data_i[:, 1]


def plot_data():
    for i in range(num_files):
        x = wavelengths
        y = dust_emission_data[i, :]
        plot.plot(x, y, linewidth = 1)

    plot.xscale('log')
    #plot.yscale('log')

    plot.show()

def interpolate_data(wavelength, absorption_fraction):
    x = wavelengths
    y = np.linspace(0, 1, num_files + 1)
    z = dust_emission_data

    dust_emission_f = interp2d(x, y, z)
    return dust_emission_f(wavelength, absorption_fraction)

read_data()
#plot_data()

test_wavelengths = np.linspace(1, 50, 10)

print interpolate_data(test_wavelengths, 0)
