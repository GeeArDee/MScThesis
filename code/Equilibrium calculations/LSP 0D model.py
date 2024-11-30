# %%
import matplotlib.pyplot as plt
import numpy as np
from CEA_LSP_polynomials import CEA_LSP_polynomials
from scipy.optimize import minimize, Bounds, fsolve
import pandas as pd
import sympy as sp
from labellines import labelLines
from Equilibrium_calculations import MinimizeGibbs_Ar, h_mix_mass_Ar


# Initialize variables
t = 0                                       # Time (s)
t_end = 10                                  # End time of the simulation (s)
n_points = 1000                             # Number of points for radiation simulation (-)

T_ini = 300.0                                # Initial temperature of the chamber (K)
p_ini = 20.0e5                                # Pressure of the chamber (20 bar to Pa)
E_laser = 30.0                                # Total input energy by the laser (J)
V_chamberV1 = 0.4/1000                      # Volume of the V1 chamber (0.4 L to m^3)                      


# Functions                                 
def rho_300K(T, p):                         # Density of argon near 300 K (kg/m^3)
    R_u = 8314.462                          # Universal gas constant [J/kmol-K]
    MW_Ar = 39.948                          # Molecular weight of argon [u]
    R_specific = R_u/MW_Ar
    return p/(R_specific*T)

def h_mix(T, p):
    # Take symbolic expression of h(T) and calculate h(T)
    t = sp.symbols('t')
    h = sp.lambdify(t, h_mix_mass_Ar(T, p), "numpy")    
    return h(T)

# def T_final(Q, m, T_1, p_1, p_2): # NO WORK, FIX THIS FUNC!
#     # Q = m * (h2 - h1), solve for T_2
#     t = sp.symbols('t')
#     h_1 = h_mix(T_1, p_1)
#     print(h_1)
#     T_guess = 6001                         # Guess that T is higher than 6000 K for the polynomials to return the right one
#     #h_2 = sp.lambdify(t, h_mix_mass_Ar(T_guess, p_2), "numpy") 
#     h_2 = h_mix_mass_Ar(T_guess, p_2)
#     print(h_2)
#     system = m*(h_2 - h_1) - Q
#     system_lambdified = sp.lambdify(t, system, "numpy")
#     #T_2 = fsolve(system, 15000)
#     return system #change to T_2

def solve_for_T2(Q, m, T_1, p_1, T_2_guess, p_2):
    def system_to_solve(T_2, T_1, p_2, p_1, Q, m):                      # Set up the equation to solve
        h_1 = h_mix(T_1, p_1)
        return m*(h_mix(T_2, p_2) - h_1) - Q
    T_2 = fsolve(system_to_solve, T_2_guess, args=(T_1,p_2,p_1,Q,m))      # Solve for T_2
    return T_2

def V_final(p, N_ions, N_electrons, T):
    # p*V = N*R_u*T
    R_u = 8314.462                          # Universal gas constant [J/kmol-K]
    N = N_ions + N_electrons
    return N*R_u*T/p


# STEP 1: Volume and mass of plasma cone at p_ini and T_ini 
l = 2/100                                   # Length of plasma cone (m)
r = 1/1000                                  # Radius of the base of plasma cone (m) CHANGE THIS, FIRST APPROX.
V_plasma = np.pi*r**2*l/3                   # Volume of plasma cone (m^3)
m_plasma = V_plasma*rho_300K(T_ini, p_ini)  # Mass of plasma (kg), will be taken as constant during this exercise


# STEP 2: Add 30 J of energy to m_plasma, while keeping constant pressure
p_2 = p_ini                                                             # Plasma cone is at the same pressure as the surrounding gas
T_2_guess = 15000.0
T_2 = solve_for_T2(E_laser, m_plasma, T_ini, p_ini, T_2_guess, p_2)     # Temperature of plasma after energy addition (K)
print(T_2)


# STEP 3: Find new volume (V_2) of cone
#V_2 = V_final(p_2, N_ions, N_electrons, T_2)


# STEP 4: Calculate pressure increase in chamber due to expansion of gas in plasma core


# STEP 5: As gas cools due to Brems. radiation, volume of the cone contracts

# %% Next steps

#while t < t_end:#

# Calculate rate of heat loss via bremsstralung.
# Is plasma ball optically thick or thin? -> volumetric/surface emitter
# Rate of heat loss? Suppose it's thick (a surface emitter)?
# delta_q = m*(h2-h1)

# %% Old code

# # Load Cp from saved file (from Variable Cp.py)
# Cp_20bar_array = np.ndarray.flatten(np.load('Cp_20bar.npy'))    # Cp of argon at 20 bar, depends on temperature (kJ/kg-K)
# T_array = np.load('T_20bar.npy')            # Temperature array going with Cp (K)
# #plt.figure()
# #plt.plot(T_array, Cp_20bar_array)
# print(T_array)
# print(Cp_20bar_array)

