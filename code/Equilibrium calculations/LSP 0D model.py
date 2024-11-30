# %%
import matplotlib.pyplot as plt
import numpy as np
from CEA_LSP_polynomials import CEA_LSP_polynomials
from scipy.optimize import minimize, Bounds, fsolve
from scipy import constants
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
MW_Ar = 39.948                          # Molecular weight of argon [u]      
Z_Ar = 18                               # Atomic number of argon [-]



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

def solve_for_T2(Q, m, T_1, p_1, T_2_guess, p_2):
    def system_to_solve(T_2, T_1, p_2, p_1, Q, m):                      # Set up the equation to solve
        h_1 = h_mix(T_1, p_1)
        return m*(h_mix(T_2, p_2) - h_1) - Q
    T_2 = fsolve(system_to_solve, T_2_guess, args=(T_1,p_2,p_1,Q,m))      # Solve for T_2
    return T_2

def V_final(p, n, T):
    # p*V = n*R_u*T
    R_u = 8314.462                          # Universal gas constant [J/kmol-K]
    return n*R_u*T/p

def P_brems_perV(T_e, n_i, n_e, Z):              # Total radiation power by Bremsstrahlung per unit volume [W/m^3 (?)] -> Check units! From Plasma physics: An introduction to laboratory, space, and fusion plasmas (2010) by A. Piel
    constant_term = 8*constants.pi**2 * (constants.k)**(1/2) / (np.sqrt(3) * (4*constants.pi*constants.epsilon_0)**3 * constants.m_e**(3/2) * constants.c**3 * constants.h)
    print(constant_term) # should be close to 5.34e-37
    return constant_term * Z**2 * n_e * n_i * T_e**(1/2)



# STEP 1: Volume and mass of plasma cone at p_ini and T_ini 
l = 2/100                                   # Length of plasma cone (m)
r = 1/1000                                  # Radius of the base of plasma cone (m) CHANGE THIS, FIRST APPROX.
V_plasma = np.pi*r**2*l/3                   # Volume of plasma cone (m^3)
m_plasma = V_plasma*rho_300K(T_ini, p_ini)  # Mass of plasma (kg), will be taken as constant during this exercise


# STEP 2: Add 30 J of energy to m_plasma, while keeping constant pressure
p_2 = p_ini                                                             # Plasma cone is at the same pressure as the surrounding gas
T_2_guess = 15000.0
T_2 = solve_for_T2(E_laser, m_plasma, T_ini, p_ini, T_2_guess, p_2)[0]     # Temperature of plasma after energy addition (K)
# print(T_2)


# STEP 3: Find new volume (V_2) of cone
n_tot_ini = m_plasma/MW_Ar                                      # Calculate initial number of moles 
(Gibbs, x, y_Ar, y_ArII, y_e) = MinimizeGibbs_Ar(T_2, p_2)         # Get dissociation x and mole fractions at temp and press of step 2
# print(MinimizeGibbs_Ar(T_2, p_2) )
n_tot = n_tot_ini*(x + 1)
n_Ar = y_Ar*n_tot
n_ArII = y_ArII*n_tot
n_e = y_e*n_tot
V_2 = V_final(p_2, n_tot, T_2)


# STEP 4: Calculate pressure increase in chamber due to expansion of gas in plasma core


# STEP 5: As gas cools due to Brems. radiation, volume of the cone contracts
P_brems = V_2 * P_brems_perV(T_2, n_ArII, n_e, Z_Ar)
# print(P_brems)

# %% Next steps

#while t < t_end:#

# Calculate rate of heat loss via bremsstralung.
# Is plasma ball optically thick or thin? -> volumetric/surface emitter
# Rate of heat loss? Suppose it's thick (a surface emitter)?
