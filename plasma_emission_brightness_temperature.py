import numpy as np
import matplotlib.pyplot as plt

T = 1.0e6                               # Coronal Temperature (K)
T_1 = 5e8                      		# Temperature of energetic electrons (K)
# n = 2.3e11                       	# Plasma density (cm^-2)
nu_p = 144e6                		# Plasma Frequency (Hz)
w = 10**(np.arange(-9.0, -2.9, 0.01))   # Normalized energy density of Langmuir waves
#w = 1e-5    				# Normalized energy density of Langmuir waves
#L_n = 7.56e11               		# Scale height of corona (cm)
L_n = 1e9                   		# Scale height of corona (cm)

# Fundamental constants
m = 9.1095e-28        # Electron mass (g)
k = 1.3803e-16          # Boltzmann constant (erg/K)
c = 2.9979e10           # Speed of light (cm/s)
mass_ratio = 1/1836.15        # Mass ratio of electrons to ions

# Calculated parameters
# nu_p = 9000 * n**(1/2)
# print("nu_p: ", "{:.2e}".format(nu_p))
n = (nu_p/9000)**2
print("n: ", "{:.2e}".format(n))
omega_p = 2 * np.pi * nu_p
print("omega_p: ", "{:.2e}".format(omega_p))
nu_ei = 5.5 * n * T**(-3/2) * np.log(1e4 * T**(3/2) * n**(-1/3)) # Electron-ion collision frequency
print("nu_ei: ", "{:.2e}".format(nu_ei))
v_T = ((k * T)/m)**(1/2)      # Thermal velocity of background electrons
print("v_T: ", "{:.2e}".format(v_T))
v_1 = c * (1 - ((m*c**2)/(k*T_1 + m*c**2))**2)**(1/2)
print("v_1: ", "{:.2e}".format(v_1))

def brightness_temp_fundamental(w, T, nu_p, L_n, n, omega_p, nu_ei, v_T, v_1):
	k_min = omega_p/c       # Minimum wave number
	k_max = omega_p/(5*v_T) # Maximum wave number

	A = (np.pi/36) * ((v_1**2)/(v_T**2)) * omega_p * T
	B = 2 * np.sqrt(3) * (v_T/c) * (k_max - k_min) * omega_p**(-1)
	C = (np.pi/108) * mass_ratio * (v_1**2/v_T**2) * omega_p
	#print("A: ", "{:.2e}".format(A))
	#print("B: ", "{:.2e}".format(B))
	#print("C: ", "{:.2e}".format(C))

	#print("{:.2e}".format(A*w))
	#print("{:.2e}".format(nu_ei - C*w))

	T_bf = ((A*w)/(nu_ei - C*w)) * (1 - np.exp(-B * (nu_ei - C*w) * L_n))
	return T_bf

brightness_temps = np.zeros(len(w))
for i in range(len(w)):
	brightness_temps[i] = brightness_temp_fundamental(w[i], T, nu_p, L_n, n, omega_p, nu_ei, v_T, v_1)

plt.figure(0)
plt.plot(np.log10(w), np.log10(brightness_temps), 'b-')
plt.xlim(-9, -3)
plt.ylim(8,18)

plt.show()

#T_bf = brightness_temp_fundamental(1e-5, T, nu_p, L_n, n, omega_p, nu_ei, v_T, v_1)
#print("T_bf: ", "{:.2e}".format(T_bf))
