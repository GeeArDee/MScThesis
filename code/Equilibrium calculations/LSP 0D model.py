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
k = 1.666                               # k, or gamma


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

def p_final(V, n, T):
    # p*V = n*R_u*T
    R_u = 8314.462                          # Universal gas constant [J/kmol-K]
    return n*R_u*T/V

def P_brems_perV(T_e_K, n_i, n_e, Z):              # Total radiation power by Bremsstrahlung per unit volume [W/m^3 (?)]
    constant_term = 1.57e-28
    return constant_term * Z**2 * n_e * n_i * T_e_K**(1/2)

def rho_N_ArII(T, p, n_tot_ini):                                # Number density of argon ions (n/m^3) = Number density of electrons (n/m^3)
    (Gibbs, x, y_Ar, y_ArII, y_e) = MinimizeGibbs_Ar(T, p)      # Get dissociation x and mole fractions at temp and press
    n_tot = n_tot_ini*(x + 1)
    n_ArII = y_ArII*n_tot
    V_2 = V_final(p, n_tot, T)
    return n_ArII * constants.Avogadro /V_2                     # Return number density of Ar+, which is equal to electron density

def area_cone(V_cone, r_fixed):
    l = 3*V_cone/(np.pi*r_fixed**2)
    A = np.pi*r_fixed*(r_fixed + np.sqrt(l**2 + r_fixed**2))
    return A



# STEP 1: Volume and mass of plasma cone at p_ini and T_ini 
l = 2/100                                   # Length of plasma cone (m)
r = 1/1000                                  # Radius of the base of plasma cone (m) CHANGE THIS, FIRST APPROX.
V_plasma = np.pi*r**2*l/3                   # Volume of plasma cone (m^3)
m_plasma = V_plasma*rho_300K(T_ini, p_ini)  # Mass of plasma (kg), will be taken as constant during this exercise


# STEP 2: Add 30 J of energy to m_plasma, while keeping constant pressure
p_2 = p_ini                                                             # Plasma cone is at the same pressure as the surrounding gas
T_2_guess = 15000.0
T_2 = solve_for_T2(E_laser, m_plasma, T_ini, p_ini, T_2_guess, p_2)[0]     # Temperature of plasma after energy addition (K)
print("T_2 = ", T_2)


# STEP 3: Find new volume (V_2) of cone
n_tot_ini = m_plasma/MW_Ar                                      # Calculate initial number of moles 
(Gibbs, x, y_Ar, y_ArII, y_e) = MinimizeGibbs_Ar(T_2, p_2)         # Get dissociation x and mole fractions at temp and press of step 2
# print(MinimizeGibbs_Ar(T_2, p_2) )
n_tot = n_tot_ini*(x + 1)
n_Ar = y_Ar*n_tot
n_ArII = y_ArII*n_tot
n_e = y_e*n_tot
V_2 = V_final(p_2, n_tot, T_2)
print("V_2 =", V_2)


# STEP 4: Calculate pressure increase in chamber due to expansion of gas in plasma core
p_4 = p_ini * ((V_chamberV1 - V_plasma)/(V_chamberV1 - V_2))**k
print("p_4 =", p_4)
print("Pressure increase :", (p_4-p_ini)/1000,"kPa")

# STEP 5: As gas cools due to Brems. radiation, volume of the cone contracts
P_brems = V_2 * P_brems_perV(T_2, rho_N_ArII(T_2, p_ini, n_tot_ini), rho_N_ArII(T_2, p_ini, n_tot_ini), 1)
print("P_brems =", "{:e}".format(P_brems), "W")

    # Graph Brems. radiation with temperature
t_vector = np.arange(500, 20001, 500, dtype=float)  # temperature in K
P_brems_vector = np.zeros(np.shape(t_vector))
for i, t in enumerate(t_vector):
    P_brems_vector[i] = V_2 * P_brems_perV(t, rho_N_ArII(t, p_ini, n_tot_ini), rho_N_ArII(t, p_ini, n_tot_ini), 1)
plt.figure()
plt.plot(t_vector, P_brems_vector)
plt.xlabel('Temperature [K]')
plt.ylabel('P_brems [W]')
plt.xlim([0,20000])
plt.ylim([0,5e11])
plt.title("P_brems evolution with temperature")

    # Calculate plasma frequency to compare to Brems. light emission
omega_p = np.sqrt(rho_N_ArII(T_2, p_ini, n_tot_ini)*constants.e**2/(constants.epsilon_0*constants.m_e))
print("The plasma frequency is", "{:e}".format(omega_p), "Hz")
    # If we consider visible light and above (>1e14 Hz), it's a volume emitter!


# STEP 6: Plot pressure with time, once the laser switches off

# Initialize variables
timestep = 1e-11
time = np.array(0)  # time [s]
p = np.array(p_4)   # pressure [Pa]
T = np.array(T_2)   # Temperature [K]
V = np.array(V_2)   # Volume of plasma cone [m^3]
i = 1               # iteration index

while i < 50: #p_4 > p_ini:
    E_laser = E_laser - P_brems * timestep          # Energy in the plasma ("E_laser") goes down as energy is radiated 
    
    # STEP 6.2: Have XX J of energy to m_plasma, while keeping constant pressure
    p_2 = p_ini                                                             # Plasma cone is at the same pressure as the surrounding gas
    T_2_guess = 15000.0
    T_2 = solve_for_T2(E_laser, m_plasma, T_ini, p_ini, T_2_guess, p_2)[0]  # Temperature of plasma after energy addition (K)

    # STEP 6.3: Volume of the cone contracts; find new volume (V_2) of cone
    n_tot_ini = m_plasma/MW_Ar                                          # Calculate initial number of moles 
    (Gibbs, x, y_Ar, y_ArII, y_e) = MinimizeGibbs_Ar(T_2, p_2)          # Get dissociation x and mole fractions at temp and press of step 2
    n_tot = n_tot_ini*(x + 1)
    V_2 = V_final(p_2, n_tot, T_2)
    
    # STEP 6.4: Calculate pressure increase in chamber due to expansion of gas in plasma core
    p_4 = p_ini * ((V_chamberV1 - V_plasma)/(V_chamberV1 - V_2))**k

    # STEP 6.5: As gas cools due to Brems. radiation, volume of the cone contracts
    P_brems = V_2 * P_brems_perV(T_2, rho_N_ArII(T_2, p_ini, n_tot_ini), rho_N_ArII(T_2, p_ini, n_tot_ini), 1)

    # Save all values to vectors and increment iterator
    time = np.append(time, i * timestep)
    p = np.append(p, p_4)
    T = np.append(T, T_2)             # Store temperature in vector for plotting
    V = np.append(V, V_2)

    i += 1

# plot 
plt.figure()
plt.plot(time, p)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
#plt.xlim([0,20000])
#plt.ylim([0,5e11])
plt.title("Pressure change over time when laser turns off, with bremsstrahlung radiation loss")

#%% STEP 7: Pressure with time, with t=0 laser turns on. This is to see what happens as pressure rises.

# Initialize variables
timestep = 10e-6       # 1 us timestep
time = np.array([0])    # time [s]
p = np.array([p_ini])       # pressure [Pa]
T = np.array([T_ini])       # Temperature [K]
V = np.array([0])           # Volume of plasma cone [m^3]
E_plasma_array = np.array([0])    # Energy in the plasma [J]
E_plasma = 0
i = 1                   # iteration index
time_end = 10e-3        # end of the sim at 10 ms, when laser pulse is done
P_loss = 0             # Brems dissipation power [W]
P_laser = 3000          # laser power (3 kW) [W]

loss = "brems" # "brems" or "blackbody"

time = np.append(time, i*timestep)  # time of first iteration

while time[i] < time_end:
    E_plasma = E_plasma + P_laser * timestep - P_loss * timestep          # Energy in the plasma ("E_plasma") goes up as time progresses
    
    # STEP 6.2: Have XX J of energy to m_plasma, while keeping constant pressure
    p_2 = p_ini                                                             # Plasma cone is at the same pressure as the surrounding gas
    T_2_guess = T[i-1]      # Take the last calculated temperature as the guess for polynomials
    T_2 = solve_for_T2(E_plasma, m_plasma, T_ini, p_ini, T_2_guess, p_2)[0]  # Temperature of plasma after energy addition (K)

    # STEP 6.3: Volume of the cone contracts; find new volume (V_2) of cone
    n_tot_ini = m_plasma/MW_Ar                                          # Calculate initial number of moles 
    (Gibbs, x, y_Ar, y_ArII, y_e) = MinimizeGibbs_Ar(T_2, p_2)          # Get dissociation x and mole fractions at temp and press of step 2
    n_tot = n_tot_ini*(x + 1)
    V_2 = V_final(p_2, n_tot, T_2)
    
    # STEP 6.4: Calculate pressure increase in chamber due to expansion of gas in plasma core
    p_4 = p_ini * ((V_chamberV1 - V_plasma)/(V_chamberV1 - V_2))**k

    match loss:
        case "brems":
            # STEP 6.5: As gas heats due to Brems. radiation, volume of the cone increases
            P_loss = V_2 * P_brems_perV(T_2, rho_N_ArII(T_2, p_ini, n_tot_ini), rho_N_ArII(T_2, p_ini, n_tot_ini), 1)
            print(P_loss)
        case "blackbody":
            # calculate blackbody radiation loss, as a sanity check upper bound on power loss
            A = area_cone(V_2, r)
            P_loss = constants.Stefan_Boltzmann * T_2**4 *A # blackbody radiation (emissivity e = 1)
            print(P_loss)

    # Save all values to vectors and increment time + iterator
    p = np.append(p, p_4)
    T = np.append(T, T_2)               # Store temperature in vector for plotting
    V = np.append(V, V_2)
    E_plasma_array = np.append(E_plasma_array, E_plasma)

    i += 1
    time = np.append(time, i*timestep)  # time of next iteration

# plot 

plt.figure()
plt.plot(time, p)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
#plt.xlim([0,20000])
#plt.ylim([0,5e11])
plt.title("Pressure change over time during laser pulse, with bremsstrahlung radiation loss")
# %%
