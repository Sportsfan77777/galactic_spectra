"""
Galaxy Spectrum Tool
"""

import math
import numpy as np
from scipy.interpolate import interp2d

from matplotlib import pyplot as plot
from matplotlib import gridspec as grid
from matplotlib.widgets import Slider, Button, RadioButtons

###############################################################################

# Constants
h = 6.626 * 10**(-34)
c = 3.0 * 10**(8)
k_B = 1.38 * 10**(-23)

###############################################################################

### GUI: Figure and Subplots ###
fig = plot.figure(figsize = (10, 8))
fig.canvas.set_window_title("Galaxy Spectra Tool")

ax = plot.subplot2grid((1, 4), (0, 0), colspan = 2)
ax2 = plot.subplot2grid((1, 4), (0, 2), colspan = 2)

ax.get_shared_y_axes().join(ax, ax2)

fig.subplots_adjust(wspace = 0.50)
#ax = axes[0]; ax2 = axes[1]
plot.subplots_adjust(bottom = 0.40)

### Plotting parameters ###
fontsize = 18
linewidth = 4

normalize_to_attenuated = [False]

### Axes ###

ax.set_xlim(100, 700) #; ax.set_xscale('log')
ax2.set_xlim(100, 1000000); ax2.set_xscale('log')
ax.set_ylim(10**-4, 1)

#ax.set_xlabel("UV    Visible Light            ", fontsize = fontsize)
#ax2.set_xlabel("Infrared (IR)", fontsize = fontsize)

ax.set_xlabel("Wavelength (nm)", fontsize = fontsize)
ax2.set_xlabel("Wavelength (nm)", fontsize = fontsize)

ax.set_ylabel("Luminosity (normalized)", fontsize = fontsize)

ax.set_title("Zoom-in\n(UV and Visible)", fontsize = fontsize + 2)
ax2.set_title("UV (<400 nm)\nVisible (400 to 700 nm)\n IR (>700 nm)", fontsize = fontsize - 3)

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

def planck(wavelength, temperature):
    # Planck's Law
    wavelength_m = wavelength * 10**(-9) # convert nm to m
    return np.power(wavelength_m, -5.0) / np.expm1(h * c / k_B / wavelength_m / temperature)

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

### Helper functions ###

###############################################################################

# Initial plot
optical_wavelengths = np.linspace(50, 700, 1000) 
ir_wavelengths = np.logspace(np.log10(700), np.log10(1000000), 10000)
wavelengths = np.concatenate((optical_wavelengths, ir_wavelengths))

fluxes = planck(wavelengths, 15000)
fluxes /= np.max(fluxes)
spectrum, = ax.plot(wavelengths, fluxes, c = "b", linewidth = linewidth)
spectrum2, = ax2.plot(wavelengths, fluxes, c = "b", linewidth = linewidth)

spectrum_extincted, = ax.plot(wavelengths, fluxes, c = "r", linewidth = linewidth)
spectrum_extincted2, = ax2.plot(wavelengths, fluxes, c = "r", linewidth = linewidth)

print wavelengths
print fluxes

#### Reference lines ####
y_ref = [10**(-6), 1.5]
ax.plot([400, 400], y_ref, c = "k", linewidth = linewidth - 1)
ax.plot([700, 700], y_ref, c = "k", linewidth = linewidth + 1)

ax2.plot([400, 400], y_ref, c = "k", linewidth = linewidth - 1)
ax2.plot([700, 700], y_ref, c = "k", linewidth = linewidth - 1)

###############################################################################

#### Sliders ####

slider_x = 0.42
slider_y = 0.25
slider_length = 0.56 - slider_x
slider_height = 0.02
slider_separation = 0.03

ax_b = plot.axes([slider_x, slider_y - 1.0 * slider_separation, slider_length, slider_height])
ax_a = plot.axes([slider_x, slider_y - 2.0 * slider_separation, slider_length, slider_height])
ax_g = plot.axes([slider_x, slider_y - 3.0 * slider_separation, slider_length, slider_height])
ax_m = plot.axes([slider_x, slider_y - 4.0 * slider_separation, slider_length, slider_height])

b_slider = Slider(ax_b, 'B stars (Blue)', -1, 5.0, valinit = 0)
a_slider = Slider(ax_a, 'A stars (White)', -1, 5.0, valinit = -0.01)
g_slider = Slider(ax_g, 'G stars (Yellow)', -1, 5.0, valinit = -0.01)
m_slider = Slider(ax_m, 'M stars (Red)', -1, 5.0, valinit = -0.01)

T_slider_x = 0.77
T_slider_length = 0.08

ax_tb = plot.axes([T_slider_x, slider_y - 1.0 * slider_separation, T_slider_length, slider_height])
ax_ta = plot.axes([T_slider_x, slider_y - 2.0 * slider_separation, T_slider_length, slider_height])
ax_tg = plot.axes([T_slider_x, slider_y - 3.0 * slider_separation, T_slider_length, slider_height])
ax_tm = plot.axes([T_slider_x, slider_y - 4.0 * slider_separation, T_slider_length, slider_height])

Tb_slider = Slider(ax_tb, 'T (K)', 10000, 30000, valinit = 15000, valfmt = "%d")
Ta_slider = Slider(ax_ta, 'T (K)', 7500, 10000, valinit = 8000, valfmt = "%d")
Tg_slider = Slider(ax_tg, 'T (K)', 5200, 6000, valinit = 5500, valfmt = "%d")
Tm_slider = Slider(ax_tm, 'T (K)', 2400, 3700, valinit = 3000, valfmt = "%d")

ax_dust = plot.axes([slider_x, slider_y - 6.0 * slider_separation, slider_length, slider_height])

dust_slider = Slider(ax_dust, 'Dust', -1, 2.5, valinit = 0)

###############################################################################

#### Radio Buttons ####

# Radio Buttons
radio_length = 0.10
radio_height = 0.10
radio_separation = 0.12
ax_y = plot.axes([0.05, slider_y + slider_height - radio_height, radio_length, radio_height])
radio_y = RadioButtons(ax_y, ('linear', 'log'), active = 0)

ax_norm = plot.axes([0.05, slider_y + slider_height - radio_height - radio_separation, radio_length + 0.10, radio_height])
radio_norm = RadioButtons(ax_norm, ('normalize blue', 'normalize red'), active = 0)

###############################################################################

### Slider Listeners ###

def update(val):
    num_b = b_slider.val; num_a = a_slider.val; num_g = g_slider.val; num_m = m_slider.val
    temp_b = Tb_slider.val; temp_a = Ta_slider.val; temp_g = Tg_slider.val; temp_m = Tm_slider.val
    dust_Av = dust_slider.val

    include_b = 1; include_a = 1; include_g = 1; include_m = 1

    if num_b < 0:
        include_b = 0
    if num_a < 0:
        include_a = 0
    if num_g < 0:
        include_g = 0
    if num_m < 0:
        include_m = 0
    if dust_Av < 0:
        dust_Av = 0

    ### Unattenuated spectrum ###
    composite_spectrum_f = lambda x : include_b * np.power(10.0, num_b) * planck(x, temp_b) \
                                    + include_a * np.power(10.0, num_a) * planck(x, temp_a) \
                                    + include_g * np.power(10.0, num_g) * planck(x, temp_g) \
                                    + include_m * np.power(10.0, num_m) * planck(x, temp_m)
    composite_spectrum = composite_spectrum_f(wavelengths)
    composite_spectrum /= np.max(composite_spectrum)

    ### Attenuated spectrum ###
    extinction_f = np.vectorize(lambda x : dust_extinction(x, dust_Av))
    composite_spectrum_extincted = composite_spectrum * extinction_f(wavelengths)

    # Add dust emission ###
    absorption_fraction = np.sum(composite_spectrum_extincted) / np.sum(composite_spectrum)
    print "Absorption Fraction: %.3f" % absorption_fraction

    dust_emission_f = np.vectorize(lambda x : dust_emission(x, absorption_fraction, dust_emission_interpolator))
    dust_spectrum = dust_emission_f(wavelengths)

    composite_spectrum_extincted += dust_spectrum

    ### Switch normalization? ###
    if normalize_to_attenuated[0]:
        max_ratio = np.max(composite_spectrum) / np.max(composite_spectrum_extincted)
        composite_spectrum *= max_ratio
        composite_spectrum_extincted /= np.max(composite_spectrum_extincted)        
    
    spectrum.set_ydata(composite_spectrum); spectrum2.set_ydata(composite_spectrum) 
    spectrum_extincted.set_ydata(composite_spectrum_extincted); spectrum_extincted2.set_ydata(composite_spectrum_extincted)

    fig.canvas.draw_idle()

b_slider.on_changed(update); a_slider.on_changed(update); g_slider.on_changed(update); m_slider.on_changed(update)
Tb_slider.on_changed(update); Ta_slider.on_changed(update); Tg_slider.on_changed(update); Tm_slider.on_changed(update)

dust_slider.on_changed(update)

###############################################################################

### Radio Listeners ###

def set_y(val):
    ax.set_yscale(val)
    ax2.set_yscale(val)
    fig.canvas.draw_idle()

def set_norm(val):
    if val == "normalize blue":
        normalize_to_attenuated[0] = False
    if val == "normalize red":
        normalize_to_attenuated[0] = True
    update(val)
    fig.canvas.draw_idle()

radio_y.on_clicked(set_y)
radio_norm.on_clicked(set_norm)

###############################################################################

### MAIN ###
plot.show()
