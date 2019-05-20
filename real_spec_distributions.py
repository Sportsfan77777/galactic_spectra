"""
Galaxy Spectrum Tool
"""

import math
import numpy as np
from scipy.interpolate import interp1d, interp2d

from matplotlib import pyplot as plot
from matplotlib import gridspec as grid
from matplotlib.widgets import Slider, Button, RadioButtons

###############################################################################

# Constants
h = 6.626 * 10**(-34)
c = 3.0 * 10**(8)
k_B = 1.38 * 10**(-23)

# Reference wavelength and temperatures
log_stellar_temperatures = {} # Taken from the paper
log_stellar_temperatures['O'] = 4.55; log_stellar_temperatures['B'] = 4.20; # B should be 4.15
log_stellar_temperatures['A'] = 3.93; log_stellar_temperatures['F'] = 3.83; log_stellar_temperatures['G'] = 3.75
log_stellar_temperatures['K'] = 3.69; log_stellar_temperatures['M'] = 3.55;

stellar_temperatures = {} # Same as above, but more readable
stellar_temperatures['O'] = 31622; stellar_temperatures['B'] = 15848; # B should be 14125
stellar_temperatures['A'] = 8511; stellar_temperatures['F'] = 6761; stellar_temperatures['G'] = 5623
stellar_temperatures['K'] = 4897; stellar_temperatures['M'] = 3548;

ref_wavelength = 555.6
ref_temperature = stellar_temperatures['M']

###############################################################################

### Initial Mass Function ###

# Using \alpha = 2.35 for M > 0.5 and \alpha = 1.35 for M < 0.5
initial_mass_function = np.array([0.016, 0.254, 0.198, 0.232, 0.299, 1.160, 6.275]) # OBAFGKM

###############################################################################

### Interpolators ###

def get_radius_interpolator():
    # Using reference values from https://en.wikipedia.org/wiki/Main_sequence
    temperatures = np.array([2660, 3120, 3920, 4410, 5240, 5610, 5780, 5920, 6540, 7240, 8620, 10800, 16400, 30000, 38000])
    radii = np.array([0.13, 0.32, 0.63, 0.74, 0.85, 0.93, 1.0, 1.05, 1.2, 1.3, 1.7, 2.5, 3.8, 7.4, 18.0])

    interpolator = interp1d(temperatures, radii)
    return interpolator
radius_interpolator = get_radius_interpolator()

def O_star_interpolator():
    data = np.loadtxt("stellar_spectra/O9_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
O_interpolator = O_star_interpolator()

def B_star_interpolator():
    data = np.loadtxt("stellar_spectra/B57_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
B_interpolator = B_star_interpolator()

def A_star_interpolator():
    data = np.loadtxt("stellar_spectra/A5_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
A_interpolator = A_star_interpolator()

def F_star_interpolator():
    data = np.loadtxt("stellar_spectra/F2_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
F_interpolator = F_star_interpolator()

def G_star_interpolator():
    data = np.loadtxt("stellar_spectra/G2_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
G_interpolator = G_star_interpolator()

def K_star_interpolator():
    data = np.loadtxt("stellar_spectra/K2_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
K_interpolator = K_star_interpolator()

def M_star_interpolator():
    data = np.loadtxt("stellar_spectra/M2_star.txt")
    wavelengths = data[:, 0]; spectrum = data[:, 1]
    return interp1d(wavelengths, spectrum)
M_interpolator = M_star_interpolator()

stellar_spectra_interpolators = {}
stellar_spectra_interpolators['O'] = O_interpolator; stellar_spectra_interpolators['B'] = B_interpolator
stellar_spectra_interpolators['A'] = A_interpolator; stellar_spectra_interpolators['F'] = F_interpolator; stellar_spectra_interpolators['G'] = G_interpolator
stellar_spectra_interpolators['K'] = K_interpolator; stellar_spectra_interpolators['M'] = M_interpolator

###############################################################################

### GUI: Figure and Subplots ###
fig, ax = plot.subplots(figsize = (10, 6))
fig.canvas.set_window_title("Galaxy Spectra Tool")

plot.subplots_adjust(bottom = 0.45)

### Plotting parameters ###
fontsize = 16
linewidth = 2

normalize_to_attenuated = [False]

### Axes ###

ax.set_xlim(100, 1000) #; ax.set_xscale('log')
ax.set_ylim(10**-4, 1)

#ax.set_xlabel("UV    Visible Light            ", fontsize = fontsize)
#ax2.set_xlabel("Infrared (IR)", fontsize = fontsize)

ax.set_xlabel("Wavelength (nm)", fontsize = fontsize)
ax.set_ylabel("Luminosity (normalized)", fontsize = fontsize)

#ax.set_title("Zoom-in\n(UV, Visible, and Near-IR)", fontsize = fontsize + 2)
#x2.set_title("UV (10 to 400 nm)\nVisible (400 to 700 nm)\n IR (700 to 10^6 nm)", fontsize = fontsize - 3)

#plot.xscale('log')

###############################################################################

### Plotting functions ###

def get_dust_interpolation_function():
    # From Dale, Helou, et al. 2014
    num_files = 64

    dust_emission_data = np.zeros((num_files + 1, 215)) # There are 214 wavelengths in each data file
    dust_emission_data[:-1, 0] = -35; dust_emission_data[-1, :] = -50 # Have a zero absorption case

    dust_wavelengths = np.zeros(215)
    dust_wavelengths[0] = 0.05

    for i in range(num_files):
        fn = "dust_emission/spec_%02d.dat" % (i + 1)
        data_i = np.loadtxt(fn)
        dust_wavelengths[1:] = data_i[:, 0]
        dust_emission_data[i, 1:] = data_i[:, 1]

    x = dust_wavelengths * 1000.0 # Convert from um to nm
    y = np.linspace(0, 1, num_files + 1) # Absorption fraction
    z = dust_emission_data - np.max(dust_emission_data) - 1.0

    dust_emission_interpolator = interp2d(x, y, z)
    return dust_emission_interpolator
dust_emission_interpolator = get_dust_interpolation_function()

### Spectrum functions ###

def planck(wavelength, stellar_type):
    # Wavelength range for interpolator
    start_wavelength = 150.0; end_wavelength = 2450.0
    if stellar_type == 'M':
        start_wavelength = 300.0
    interpolator = stellar_spectra_interpolators[stellar_type]

    # Temperature and Radius
    temperature = np.power(10, log_stellar_temperatures[stellar_type])
    radius = radius_interpolator(temperature)

    # Helper for normalization
    def plancks_law(wavelength, temperature):
        wavelength_meters = wavelength * 10**-9
        return np.power(wavelength_meters, -5.0) / np.expm1(h * c / k_B / wavelength_meters / temperature)
    normalization = plancks_law(ref_wavelength, temperature) / plancks_law(ref_wavelength, ref_temperature)

    # Get spectral radiance
    if start_wavelength <= wavelength and wavelength <= end_wavelength:
        # from reference spectra
        planck_value = interpolator(wavelength) * normalization
    elif wavelength < start_wavelength:
        # from Planck's law
        normalization *= interpolator(start_wavelength) / plancks_law(start_wavelength, temperature)
        planck_value = plancks_law(wavelength, temperature) * normalization
    elif wavelength > end_wavelength:
        # from Planck's law
        normalization *= interpolator(end_wavelength) / plancks_law(end_wavelength, temperature)
        planck_value = plancks_law(wavelength, temperature) * normalization
    return planck_value

def luminosity(wavelength, stellar_type):
    # Luminosity: L ~ R^2 * B(\lambda, temperature)
    temperature = np.power(10, log_stellar_temperatures[stellar_type])
    radius = radius_interpolator(temperature)
    return 4.0 * np.pi * np.power(radius, 2) * planck(wavelength, stellar_type)
vectorized_luminosity = np.vectorize(luminosity)

def dust_extinction(wavelength, A_v, R_v = 4.05):
    # Calzetti Extinction Curve
    # Source: http://webast.ast.obs-mip.fr/hyperz/hyperz_manual1/node10.html
    # Source: https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_4.html
    wavelength_u = wavelength * 10**(-3) # convert nm to um
    k_lambda = 0 # extinction curve
    if wavelength_u < 0.63:
        k_lambda = 2.659 * (-2.156 + 1.509 / wavelength_u - 0.198 / wavelength_u**2 + 0.011 / wavelength_u**3) + R_v
    else:
        k_lambda = 2.659 * (-1.857 + 1.040 / wavelength_u) + R_v
    return np.power(10, -0.4 * k_lambda * A_v / R_v)

def dust_emission(wavelength, absorption_fraction, dust_emission_f):
    # Use interpolator
    return np.power(10, dust_emission_f(wavelength, absorption_fraction))

###############################################################################

# Rainbow Region
alpha_rainbow = 0.10; num_colors = 500

coordinates = np.linspace(400, 700, num_colors); y_region = np.array([10**(-6), 1])
visible_spectrum = np.zeros((num_colors, 2))
visible_spectrum[:, 0] = coordinates; visible_spectrum[:, 1] = coordinates
ax.pcolormesh(coordinates, y_region, np.transpose(visible_spectrum), cmap = 'nipy_spectral', alpha = alpha_rainbow)
#ax2.pcolormesh(coordinates, y_region, np.transpose(visible_spectrum), cmap = 'nipy_spectral', alpha = 0.01)

### Initial plot ###
optical_wavelengths = np.linspace(5, 1000, 2000) 
ir_wavelengths = np.logspace(np.log10(1000), np.log10(1000000), 10000)
wavelengths = np.concatenate((optical_wavelengths, ir_wavelengths))

#fluxes = vectorized_luminosity(wavelengths, 'G')
#fluxes /= np.max(fluxes)
fluxes = np.zeros(len(wavelengths))
spectrum, = ax.plot(wavelengths, fluxes, c = "b", linewidth = linewidth)
#spectrum2, = ax2.plot(wavelengths, fluxes, c = "b", linewidth = linewidth)

spectrum_extincted, = ax.plot(wavelengths, fluxes, c = "r", linewidth = linewidth)
#pectrum_extincted2, = ax2.plot(wavelengths, fluxes, c = "r", linewidth = linewidth)

#print wavelengths
#print fluxes

###############################################################################

#### Sliders ####

slider_x = 0.42
slider_y = 0.30
slider_length = 0.56 - slider_x
slider_height = 0.02
slider_separation = 0.03

ax_pop1 = plot.axes([slider_x, slider_y - 1.0 * slider_separation, slider_length, slider_height])
ax_pop2 = plot.axes([slider_x, slider_y - 2.0 * slider_separation, slider_length, slider_height])
ax_pop3 = plot.axes([slider_x, slider_y - 3.0 * slider_separation, slider_length, slider_height])

pop1_slider = Slider(ax_pop1, 'Population #1', 0, 10, valinit = 0.0)
pop2_slider = Slider(ax_pop2, 'Population #2', 0, 10, valinit = 0.0)
pop3_slider = Slider(ax_pop3, 'Population #3', 0, 10, valinit = 0.0)

T_slider_x = 0.82
T_slider_length = 0.08

ax_T_pop1 = plot.axes([T_slider_x, slider_y - 1.0 * slider_separation, T_slider_length, slider_height])
ax_T_pop2 = plot.axes([T_slider_x, slider_y - 2.0 * slider_separation, T_slider_length, slider_height])
ax_T_pop3 = plot.axes([T_slider_x, slider_y - 3.0 * slider_separation, T_slider_length, slider_height])

T_pop1_slider = Slider(ax_T_pop1, 'Hottest star in Pop #1', 0, 7, valinit = 7, valstep = 1, valfmt = "%d")
T_pop2_slider = Slider(ax_T_pop2, 'Hottest star in Pop #2', 0, 7, valinit = 5, valstep = 1, valfmt = "%d")
T_pop3_slider = Slider(ax_T_pop3, 'Hottest star in Pop #3', 0, 7, valinit = 3, valstep = 1, valfmt = "%d")

###############################################################################

### Luminosities --- Figure out somewhere else to put this ###

luminosity_O = vectorized_luminosity(wavelengths, 'O')
luminosity_B = vectorized_luminosity(wavelengths, 'B')
luminosity_A = vectorized_luminosity(wavelengths, 'A')
luminosity_F = vectorized_luminosity(wavelengths, 'F')
luminosity_G = vectorized_luminosity(wavelengths, 'G')
luminosity_K = vectorized_luminosity(wavelengths, 'K')
luminosity_M = vectorized_luminosity(wavelengths, 'M')

#luminosity_array = [luminosity_O, luminosity_B, luminosity_A, luminosity_F, luminosity_G, luminosity_K, luminosity_M]
luminosity_array = [luminosity_M, luminosity_K, luminosity_G, luminosity_F, luminosity_A, luminosity_B, luminosity_O]

def sum_luminosities(include):
    luminosities = np.zeros(len(wavelengths))
    for (include_i, luminosity_i) in zip(include, luminosity_array):
        luminosities += (include_i * luminosity_i)
    return luminosities

###############################################################################

### Slider Listeners ###

def update(val):
    num_1 = pop1_slider.val; num_2 = pop2_slider.val; num_3 = pop3_slider.val
    temp_1 = int(T_pop1_slider.val); temp_2 = int(T_pop2_slider.val); temp_3 = int(T_pop3_slider.val)
    #dust_Av = dust_slider.val
    dust_Av = 0

    include1 = np.zeros(7); include2 = np.zeros(7); include3 = np.zeros(7) # There are 7 stellar types
    include1[:temp_1] = 1; include2[:temp_2] = 1; include3[:temp_3] = 1

    ### Unattenuated spectrum ### 
    composite_spectrum1 = (10**num_1 - 1) * sum_luminosities(include1)
    composite_spectrum2 = (10**num_2 - 1) * sum_luminosities(include2)
    composite_spectrum3 = (10**num_3 - 1) * sum_luminosities(include3)
    composite_spectrum = composite_spectrum1 + composite_spectrum2 + composite_spectrum3

    composite_spectrum /= np.max(composite_spectrum)

    ### Attenuated spectrum ###
    extinction_f = np.vectorize(lambda x : dust_extinction(x, dust_Av))
    composite_spectrum_extincted = composite_spectrum * extinction_f(wavelengths)

    # Add dust emission ###
    absorption_fraction = np.sum(composite_spectrum_extincted) / np.sum(composite_spectrum)
    #print "Absorption Fraction: %.3f" % absorption_fraction

    dust_emission_f = np.vectorize(lambda x : dust_emission(x, absorption_fraction, dust_emission_interpolator))
    dust_spectrum = dust_emission_f(wavelengths)

    composite_spectrum_extincted += dust_spectrum

    ### Switch normalization? ###
    if normalize_to_attenuated[0]:
        max_ratio = np.max(composite_spectrum) / np.max(composite_spectrum_extincted)
        composite_spectrum *= max_ratio
        composite_spectrum_extincted /= np.max(composite_spectrum_extincted)        
    
    spectrum.set_ydata(composite_spectrum)
    spectrum_extincted.set_ydata(composite_spectrum_extincted)

    fig.canvas.draw_idle()

pop1_slider.on_changed(update); pop2_slider.on_changed(update); pop3_slider.on_changed(update)
T_pop1_slider.on_changed(update); T_pop2_slider.on_changed(update); T_pop3_slider.on_changed(update)

### MAIN ###
plot.show()
