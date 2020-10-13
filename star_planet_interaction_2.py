import numpy as np
from scipy.optimize import fsolve

# Physical constants (in SI)
m_p = 1.672e-27
k_B = 1.38e-23
G = 6.674e-11
u_0 = 1.2566e-6

# Units (in SI)
solar_mass = 1.989e30		# Solar mass in kg
solar_radius = 6.957e8		# Solar radius in m
jupiter_radius = 7.1492e7	# Jupiter radius in m
au = 1.496e11                   # Astronomical unit in m
parsec = 3.086e16               # Parsec in m
jansky = 1e-26

# Set physical parameters of star and planet
T_wind = 12.6e6			# Temperature of stellar wind in K (same as coronal temp.)
M_star = 0.100 * solar_mass 	# Mass of star in solar mass units
R_star = 1.000 * solar_radius	# Radius of star in solar radii
mass_loss_rate = 1e-9 * solar_mass      # Mass loss rate of star in solar mass units per year
R_p = 1.00* jupiter_radius	# Radius of planet in Jupiter radii
r = 10 *R_star    		# Radial distance of planet from star
B_star_0 = 0.1                  # Magnetic field of star in T
B_p = 0.009        		# Max. surface mag. field strength of planet
P_s = 1.86      		# Rotation period of star 
d = 151.9*parsec		# Distance to star in pc

# Set parameters of observed radio emission
flux_density = 28.15e-3 * jansky          # Flux density in Jy
solid_angle = 1.6                   # Beam solid angle in steradians
rad_efficiency = 0.01

bandwidth_MHz = 2.8 * 1e4 * B_star_0
print("Bandwidth (MHz): ", bandwidth_MHz)

# Convert period to seconds
P_s = P_s * 24 * 60 * 60
# Convert mass loss rate to per second
mass_loss_rate = mass_loss_rate / (365*24*60*60)
# Convert bandwidth to Hz
bandwidth = bandwidth_MHz*1e6       

c_s = np.sqrt((k_B* T_wind)/m_p)
r_c = (m_p*G*M_star) / (4*k_B*T_wind)

print("c_s: ", c_s)
print("r_c: ", r_c)

def parker_solar_wind(v_wind):
	return (v_wind**2/c_s**2) - np.log(v_wind**2/c_s**2) - 4*np.log((r-R_star)/r_c) - 4*(r_c/(r-R_star)) + 3

v_wind = fsolve(parker_solar_wind, c_s)[0]
print("v_wind: ", v_wind)

# Keplerian velocity
v_k = np.sqrt(M_star*G/r)
print("v_k: ", v_k)

# Effective velocity (Stellar wind velocity from planet frame of rest)
v_eff = np.sqrt(v_k**2 + v_wind**2)
print("v_eff: ", v_eff)

# Stellar wind density
density_sw = mass_loss_rate / (4*np.pi*v_wind*(r**2))
print("Mass density (kg/m^3): ", density_sw)

Omega_s = (2*np.pi) / P_s			# Ang. vel. of star
B_r = B_star_0 * (R_star / r)**2	# Radial mag. field
B_phi = B_r * ((Omega_s*r)/v_eff)		# Azim. mag. field

print("B_r (G): ", B_r*1e4)
print("B_phi (G): ", B_phi*1e4)

theta = np.arctan(B_phi/B_r) - np.arctan(v_k/v_wind)
B_star = np.sqrt(B_r**2 + B_phi**2)
B_perp = B_star * np.sin(theta)

print("B_star (G): ", B_star*1e4)
print("B_perp (G): ", B_perp*1e4)

# Alfven speed
v_A = B_r / np.sqrt(u_0 * density_sw)
print("Alfven speed: (m/s): ", v_A)
# Alfven Mach number
M_A = v_eff / v_A
print("Alfven Mach number: ", M_A)

# Calculate radius of magnetosphere
R_m = ( ( 2.44**2 * B_p**2 ) / (2*u_0 * ( density_sw*v_eff**2 + (density_sw*k_B*T_wind/m_p) + ((B_star**2)/(2*u_0)) ) ) )**(1/6) * R_p

print("R_m/R_p: ", R_m/R_p)
if R_m/R_p < 1.0:
    R_m = R_p

# Electical field of stellar wind
E_sw = v_wind * B_star * np.sin( np.arctan(B_phi/B_r) )
print("Electrical field of stellar wind (V/m): ", E_sw)

################################################################################
############ Calculate theoretical Poynting flux density on star ###############
################################################################################
#S_th = 2*np.pi*R_m**2 * (E_sw*B_perp/u_0) * M_A
S_th = 2*np.pi*R_m**2 * ((v_eff * (B_star * np.sin(theta))**2)/u_0) * M_A
print("Theoretical Poynting flux density (W): ", S_th)
S_th_max = 2*np.pi*R_m**2 * (E_sw*B_star/u_0) * M_A
print("Theoretical max. Poynting flux density (W): ", S_th_max)

################################################################################
### Calculate observationally inferred Poynting flux density on star ###########
################################################################################
total_radio_power = flux_density * solid_angle * d**2 * bandwidth
print("Total emitted radio power (W): ", total_radio_power)

obs_poynting_flux = total_radio_power/rad_efficiency
print("Observationally inferred star-ward Poynting flux (W): ", obs_poynting_flux)

