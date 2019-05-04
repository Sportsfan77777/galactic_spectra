"""
Galaxy Spectrum Tool
"""

import math
import numpy as np
from matplotlib import pyplot as plot
from matplotlib.widgets import Slider, Button, RadioButtons

# GUI
fig, axes = plot.subplots(1, 2, sharey = True, figsize = (8, 8))
ax = axes[0]; ax2 = axes[1]
plot.subplots_adjust(bottom = 0.40)

# Constants
h = 6.626 * 10**(-34)
c = 3.0 * 10**(8)
k_B = 1.38 * 10**(-23)

# Plotting parameters
fontsize = 16
linewidth = 4

normalize_to_attenuated = [False]

# Plotting

# Spectrum functions
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

def dust_emission(wavelength):
    #Source: https://www.astro.umd.edu/~richard/ASTRO421/Lecture_9.pdf
    pass

# Initial plot
optical_wavelengths = np.linspace(50, 700, 1000) 
ir_wavelengths = np.logspace(np.log10(700), np.log10(10000), 1000)
wavelengths = np.concatenate((optical_wavelengths, ir_wavelengths))

fluxes = planck(wavelengths, 15000)
fluxes /= np.max(fluxes)
spectrum, = plot.plot(wavelengths, fluxes, linewidth = linewidth)

spectrum_extincted, = plot.plot(wavelengths, fluxes, c = "r", linewidth = linewidth)

print wavelengths
print fluxes

# Reference lines
y_ref = [10**(-6), 1.5]
plot.plot([400, 400], y_ref, c = "k", linewidth = linewidth - 1)
plot.plot([700, 700], y_ref, c = "k", linewidth = linewidth - 1)

plot.xlim(100, 1000)
plot.ylim(10**-4, 1)

plot.xlabel("Wavelength (nm)", fontsize = fontsize)
plot.ylabel("Normalized Luminosity", fontsize = fontsize)

#plot.xscale('log')

# Sliders

slider_x = 0.42
slider_y = 0.25
slider_length = 0.60 - slider_x
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

T_slider_x = 0.75
T_slider_length = 0.10

ax_tb = plot.axes([T_slider_x, slider_y - 1.0 * slider_separation, T_slider_length, slider_height])
ax_ta = plot.axes([T_slider_x, slider_y - 2.0 * slider_separation, T_slider_length, slider_height])
ax_tg = plot.axes([T_slider_x, slider_y - 3.0 * slider_separation, T_slider_length, slider_height])
ax_tm = plot.axes([T_slider_x, slider_y - 4.0 * slider_separation, T_slider_length, slider_height])

Tb_slider = Slider(ax_tb, 'T (K)', 10000, 30000, valinit = 15000)
Ta_slider = Slider(ax_ta, 'T (K)', 7500, 10000, valinit = 8000)
Tg_slider = Slider(ax_tg, 'T (K)', 5200, 6000, valinit = 5500)
Tm_slider = Slider(ax_tm, 'T (K)', 2400, 3700, valinit = 3000)

ax_dust = plot.axes([slider_x, slider_y - 6.0 * slider_separation, slider_length, slider_height])

dust_slider = Slider(ax_dust, 'Dust', -1, 2.0, valinit = 0)

# Radio Buttons
radio_length = 0.10
radio_height = 0.10
radio_separation = 0.12
ax_y = plot.axes([0.05, slider_y + slider_height - radio_height, radio_length, radio_height])
radio_y = RadioButtons(ax_y, ('linear', 'log'), active = 0)

ax_norm = plot.axes([0.05, slider_y + slider_height - radio_height - radio_separation, radio_length + 0.10, radio_height])
radio_norm = RadioButtons(ax_norm, ('normalize blue', 'normalize red'), active = 0)


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

    # Unattenuated spectrum
    composite_spectrum_f = lambda x : include_b * np.power(10.0, num_b) * planck(x, temp_b) \
                                    + include_a * np.power(10.0, num_a) * planck(x, temp_a) \
                                    + include_g * np.power(10.0, num_g) * planck(x, temp_g) \
                                    + include_m * np.power(10.0, num_m) * planck(x, temp_m)
    composite_spectrum = composite_spectrum_f(wavelengths)
    composite_spectrum /= np.max(composite_spectrum)

    # Attenuated spectrum
    extinction_f = np.vectorize(lambda x : dust_extinction(x, dust_Av))
    composite_spectrum_extincted = composite_spectrum * extinction_f(wavelengths)

    # Switch normalization?
    if normalize_to_attenuated[0]:
        # Doesn't work?????
        max_ratio = np.max(composite_spectrum) / np.max(composite_spectrum_extincted)
        composite_spectrum *= max_ratio
        composite_spectrum_extincted /= np.max(composite_spectrum_extincted)        
    
    spectrum.set_ydata(composite_spectrum)
    spectrum_extincted.set_ydata(composite_spectrum_extincted)

    fig.canvas.draw_idle()

b_slider.on_changed(update); a_slider.on_changed(update); g_slider.on_changed(update); m_slider.on_changed(update)
Tb_slider.on_changed(update); Ta_slider.on_changed(update); Tg_slider.on_changed(update); Tm_slider.on_changed(update)

dust_slider.on_changed(update)

def set_y(val):
    ax.set_yscale(val)
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

### MAIN ###
#makeGUI()

plot.show()