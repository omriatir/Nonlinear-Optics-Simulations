from eddington import fitting_function
import numpy as np


@fitting_function(n=3)
def squared_sinc(a, x):
    x = np.radians(x)  # The func' receives the angles in degrees and the numpy trig' func's work with radians.
    # Input Data:
    no1w = 1.6614  # no(794.538786nm) = 1.6615
    no2w = 1.6934  # no(397.281619nm) = 1.6940
    # ne_bar1w = 1.5463  # ne(794.538786nm) = 1.5463
    ne_bar2w = 1.5687  # ne(397.281619nm) = 1.5692
    ff = 794.419*10**-9  # Fundamental Wavelength (Incident beam wavelength) = 794.538786 nm
    t0 = np.radians(29.01)  # theta_0 in radians.
    I_max = 132584.93445623052  # Max value of the intensity (before normalizing).

    # Defining elements of the returned func' for simplicity.
    # a[2] is a free parameter
    # a[1] is a constant shift parameter for the impact angle
    # a[0] represents the length of the crystal L_0 (in nm according to the value of the ff wavelength)
    # We expect L_0 to be around 200,000 nm, 0.2 mm, 0.0002 m.

    snell = np.sin(x + a[1]) / no1w  # sin(alpha_i - a[1])/no(w)
    alpha_p = np.arcsin(snell)
    t_eff = t0 + alpha_p
    L_eff = a[0] / np.sqrt(1 - snell ** 2)
    ne2wt = (((np.sin(t_eff) ** 2) / ne_bar2w ** 2) + ((np.cos(t_eff) ** 2) / no2w ** 2)) ** -0.5
    sinc_square = np.sinc((2 * L_eff / ff ) * (no1w - ne2wt)) ** 2
    I_alpha = sinc_square

    return a[2] + I_alpha



