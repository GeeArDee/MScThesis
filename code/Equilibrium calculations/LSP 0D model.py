# %%
import matplotlib.pyplot as plt
import numpy as np
from CEA_LSP_polynomials import CEA_LSP_polynomials
from scipy.optimize import minimize, Bounds
import pandas as pd
import sympy as sp
from labellines import labelLines
from tabulate import tabulate


# Initialize variables
t = 0                                       # time (s)
delta_t = 1/1000                            # timestep of the simulation (s)
t_end = 10                                  # end time of simulation (s)
t_lsr_on = 0                                # time when laser is turned on (s)
t_lsr_off = 10/1000                         # time when laser is shut off (s)

T_plasma = 15000                            # Temperature of the plasma (K)
T_ini = 300                                 # Initial temperature of the chamber (K)
P_test = 20e5                               # Pressure of the chamber (20 bar to Pa)
lsr_power = 3000                            # Laser pulse power (W)
n_absorp_lsr = 0.8                          # Efficiency of laser absorption by the plasma cone (-)
delta_q = lsr_power*(t_lsr_off-t_lsr_on)    # Total input energy by the laser (s)
V_chamberV1 = 0.4/1000                      # Volume of the V1 chamber (0.4 L to m^3)

# Load Cp from saved file (from Variable Cp.py)
Cp_20bar_array = np.ndarray.flatten(np.load('Cp_20bar.npy'))    # Cp of argon at 20 bar, depends on temperature (kJ/kg-K)
T_array = np.load('T_20bar.npy')            # Temperature array going with Cp (K)
#plt.figure()
#plt.plot(T_array, Cp_20bar_array)
print(T_array)
print(Cp_20bar_array)

# Knowing cp from Equilibrium, find mass and volume of plasma ball, assuming constant pressure

# Volume of plasma cone
l = 2/100   # length of plasma cone (m)
r = 1/1000  # radius of the base of plasma cone (m) CHANGE THIS, FIRST APPROX.
V_plasma = np.pi*r**2*l/3#*(t/t_lsr_off)  # Volume of the plasma cone (m^3) grows linearly with time, reaching max volume at t_lsr_off

# %% HW problem state 3 with variable c_p
def p(m, R_specific, T, V):
    #  p V = m R_s T
    return m*R_specific*T/V

def T(Q, m, c, T_0):
    #  Q = m c delta_T
    return Q/(m*c)+T_0

V_0 = 0.0004  # m^3
MW_Ar = 39.95  # kg/kmol
p_0 = 20e5  # Pa
T_0 = 300  # K
R_u = 8314.462  # J/K-kmol
R_Ar = R_u/MW_Ar
cv_300K = 3/2 * R_Ar
cp_300K = 5/2 * R_Ar
k_300K = cp_300K/cv_300K
Q = 0.8*30  # J
V_A0 = 1/1000*V_0

# STATE 0 (Initial)
V_B0 = V_0-V_A0
m_A = p_0*V_A0/(R_Ar*T_0)
m_tot = p_0*V_0/(R_Ar*T_0)
m_B = m_tot-m_A

# STATE 3 (Heat addition via constant pressure)
# Iterate to find right variable cp
T_ini = 15000  # First guess
delta = 10
delta_max = 1  # Maximum allowable delta

while delta > delta_max:
    cp_var = np.interp(T_ini, T_array, Cp_20bar_array)
    print(cp_var)
    T_A3 = T(Q, m_A, cp_var, T_0)
    delta = abs(T_ini - T_A3)
    T_ini = T_A3

p_A3 = p_0
V_A3 = m_A*R_Ar*T_A3/p_A3
V_B3 = V_0 - V_A3
p_B3 = p_0*V_B0**k_300K/(V_B3**k_300K)
T_B3 = p_B3*V_B3/(m_B*R_Ar)

#  STATE 4 (Uniform heat addition to entire volume)
T_4 = T(Q, m_tot, cv_300K, T_0)
p_4 = p(m_tot, R_Ar, T_4, V_0)

#  Printing
table = [['STATE', 'p [kPa]', 'T [K]'],
            ['A3', p_A3/1000, T_A3],
            ['B3', p_B3/1000, T_B3],
            ['4', p_4/1000, T_4]]
print(V_A0)
print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))

# %% Next steps

#while t < t_end:

# Calculate rate of heat loss via bremsstralung.
# Is plasma ball optically thick or thin? -> volumetric/surface emitter
# Rate of heat loss? Suppose it's thick (a surface emitter)?
# delta_q = m*(h2-h1)


