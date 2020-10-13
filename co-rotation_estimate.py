import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

G = 6.674e-11                   # G constant in SI units

rot_period = 1.86               # Rotation period in days
radius = 1.0                    # Radius in solar radii
mass = 0.1                      # Mass in solar masses
B = 1000 * 1e-4                 # Magnetic field strength in T

radius = radius*6.96e8          # Convert to m
mass = mass*1.989e30            # Convert to kg

#corotation_distance = 10.0     # Distance in stellar radii at which co-rotation breakdown is occuring
#print("Co-rotation distance (m): ", corotation_distance*radius)

# Calculate angular velocity from rotation period of star
ang_vel = (2*np.pi)/(rot_period*24.*60.*60.)
spectral_lum =  1.0e10          # Spectral luminosity of flare in W
bandwidth = 48e6                # Bandwidth of emission (estimate as bandwidth of observation)

# Calculate co-rotation radius
r_c = ((G * mass) / ang_vel**2)**(1/3)
print("Co-rotation radius (m): ", r_c)
print("Co-rotation radius (R_*): ", r_c/radius)

# Calculate mag field strength at co-rotation radius
B_c = B * (radius/r_c)**3
print("Magnetic field strength at r_c (G): ", B_c*1e4)

# Calculate co-rotation velocity
vel = r_c * ang_vel
print("Velocity (m/s): ", vel)

power = spectral_lum * bandwidth # Power dissipated in flare
print("Power disspated (W): ", power)

mass_loss_rate = power / (0.5 * vel**2)
print("Mass-loss rate (kg/s) ", mass_loss_rate)

# Convert mass-loss rate to solar masses per year
mass_loss_rate = (mass_loss_rate*365*24*60*60) / 1.989e30
print("Mass-loss rate (M_sol/yr) ", mass_loss_rate)

