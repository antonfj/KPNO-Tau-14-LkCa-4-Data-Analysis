import numpy as np
from scipy.optimize import fsolve

# Physical constants (in cgs)
m_p = 1.672e-24
k_B = 1.38e-16
G = 6.674e-8

# Units (in cgs)
solar_mass = 1.989e33		# Solar mass in g
solar_radius = 6.957e10		# Solar radius in cm
jupiter_radius = 7.1492e9	# Jupiter radius in cm
au = 1.496e13                   # Astronomical unit in cm
parsec = 3.086e18               # Parsec in cm
jansky = 1e-23                  # Jansky in erg s^-1 cm^-2 Hz^-1
day = 24 * 60 * 60              # 1 day in s
year = 365 * day                # 1 year in s

# Set physical parameters of star and planet
T_wind = 12.6e6			# Temperature of stellar wind in K (same as coronal temp.)
M_star = 0.100 * solar_mass 	# Mass of star in solar mass units
R_star = 1.000 * solar_radius	# Radius of star in solar radii
mass_loss_rate = 1e-9 * (solar_mass/year)      # Mass loss rate of star in solar mass units per year
R_p = 1.00 * jupiter_radius	# Radius of planet in Jupiter radii
r = 10 * R_star    		# Radial distance of planet from star
B_star_0 = 1000                 # Magnetic field of star in G
B_p = 40          		# Max. surface mag. field strength of planet in G
P_s = 1.86 * day      		# Rotation period of star in days
d = 151.9 * parsec		# Distance to star in pc

# Set parameters of observed radio emission
flux_density = 28.15e-3 * jansky          # Flux density in Jy
solid_angle = 1.6                   # Beam solid angle in steradians
rad_efficiency = 0.01               # Radiative efficiency

bandwidth = 2.8e6 * B_star_0
print("Bandwidth (Hz): ", bandwidth)

c_s = np.sqrt((k_B* T_wind)/m_p)
r_c = (m_p*G*M_star) / (4*k_B*T_wind)

print("c_s (cm): ", c_s)
print("r_c (cm): ", r_c)

# Radial stellar wind velocity
def parker_solar_wind(v_r):
	return (v_r**2/c_s**2) - np.log(v_r**2/c_s**2) - 4*np.log((r-R_star)/r_c) - 4*(r_c/(r-R_star)) + 3

v_r = fsolve(parker_solar_wind, c_s)[0]
print("v_r (cm/s): ", v_r)

# Azimuthal stellar wind velocity
Omega_s = (2*np.pi) / P_s			# Ang. vel. of star
v_phi = (Omega_s * R_star**2)/r
print("v_phi (cm/s): ", v_phi)

# Keplerian velocity of planet
v_k = np.sqrt(M_star*G/r)
print("v_k (cm/s): ", v_k)

# Effective velocity (Stellar wind velocity from planet frame of rest)
v_eff = np.sqrt(v_k**2 + v_r**2)
print("v_eff (cm/s): ", v_eff)

# Stellar wind density
density_sw = mass_loss_rate / (4*np.pi*v_r*(r**2))
print("Mass density (g/cm^3): ", density_sw)

B_r = B_star_0 * (R_star / r)**2	# Radial mag. field
B_phi = B_r * ((v_phi - Omega_s*r)/v_r)		# Azim. mag. field

print("B_r (G): ", B_r*1e4)
print("B_phi (G): ", B_phi*1e4)

theta = np.arctan(B_phi/B_r) - np.arctan(v_k/v_r)
B_star = np.sqrt(B_r**2 + B_phi**2)
B_perp = B_star * np.sin(theta)

print("B_star (G): ", B_star*1e4)
print("B_perp (G): ", B_perp*1e4)

# Alfven speed
v_A = B_r / np.sqrt(4*np.pi * density_sw)
print("Alfven speed: (cm/s): ", v_A)
# Alfven Mach number
M_A = v_r / v_A
print("Alfven Mach number: ", M_A)

# Calculate radius of magnetosphere
R_m = ( ( 2.44**2 * B_p**2 ) / (8*np.pi * ( density_sw*v_eff**2 + (density_sw*k_B*T_wind/m_p) + ((B_star**2)/(8*np.pi)) ) ) )**(1/6) * R_p

print("R_m/R_p: ", R_m/R_p)
if R_m/R_p < 1.0:
    R_m = R_p

################################################################################
############ Calculate theoretical Poynting flux density on star ###############
################################################################################
S_th = (1/2) * R_m**2 * v_eff * (B_star * np.sin(theta))**2 * M_A
print("Theoretical Poynting flux density (erg/s): ", S_th)

################################################################################
### Calculate observationally inferred Poynting flux density on star ###########
################################################################################
total_radio_power = flux_density * solid_angle * d**2 * bandwidth
print("Total emitted radio power (erg/s): ", total_radio_power)

obs_poynting_flux = total_radio_power/rad_efficiency
print("Observationally inferred star-ward Poynting flux (erg/s): ", obs_poynting_flux)

