import numpy as np
import matplotlib.pyplot as plt

## Constants
Fundamental_Wavelength = 794.5632 * 10**-9
c = 299792458
w1 = 2*np.pi*c/Fundamental_Wavelength #angular freq
n_ordinary_w = 1.6614
n_ordinary_2w = 1.6934
n_extraordinary_w = 1.5462
n_extraordinary_2w = 1.5687
## For alpha = 0 ISH is maximal => Delta K =:0 => theta 0 in radians:
theta_0 = np.arcsin(np.sqrt((1/n_ordinary_w**2 - 1/n_ordinary_2w**2)/ (1/n_extraordinary_2w**2 - 1/n_ordinary_2w**2)))
#This is the angle between the optical axis and alpha=0 direction.

def ISH(alpha,L0):
    Angle_after_refraction = Find_Angle_Of_Refraction(n_ordinary_w,alpha)
    theta =  - Angle_after_refraction + theta_0
    Effective_L = L0/np.cos(Angle_after_refraction)
    K1 = w1*n_ordinary_w/c
    K2 = 2*w1*n_extraordinary(theta)/c
    Delta_K = 2*K1 - K2 #Phase Mismatch
    return np.sinc(Effective_L*Delta_K/np.pi/2)**2

def Find_Angle_Of_Refraction(n2,alpha): #snell's low
    return np.arcsin(np.sin(alpha)/n2)

def n_extraordinary(theta): # Calculate n_e(angle between optial axis and k vector)
    return ((np.sin(theta)**2)/n_extraordinary_2w**2+(np.cos(theta)**2)/n_ordinary_2w**2)**(-0.5)

L0_ar = np.array([0.00001,0.0001,0.001,0.01,0.1]) #our selected L0 Values in m units
alpha_ar_deg = np.linspace(-4, 4, 1000)
alpha_ar_rad = np.radians(alpha_ar_deg)


L0=L0_ar[0]
plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=0.01[mm]")
L0=L0_ar[1]
plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=0.1[mm]")
L0=L0_ar[2]
plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=1[mm]")
L0=L0_ar[3]
plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=10[mm]")
L0=L0_ar[4]
plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=100[mm]")

plt.legend()
plt.xlabel("Incident Angle Alpha[deg]")
plt.ylabel("ISH Intensity Normalized")
plt.grid()
plt.show()

#more close L0 Values:

L0 = 25*10**-5
plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=0.0568[mm]")
#L0=0.0002
#plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=0.2[mm]")
#L0=0.0003
#plt.plot(alpha_ar_deg,ISH(alpha_ar_rad,L0),label = "L0=0.3[mm]")
plt.legend()
plt.xlabel("Incident Angle Alpha[deg]")
plt.ylabel("ISH Intensity Normalized")
plt.grid()
plt.show()