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
import time

#%%
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

#%%
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


# STEP 4: Calculate pressure increase in chamber due to expansion of gas in plasma cone
p_4 = p_ini * ((V_chamberV1 - V_plasma)/(V_chamberV1 - V_2))**k
print("p_4 =", p_4)
print("Pressure increase :", (p_4-p_ini)/1000,"kPa")

# STEP 5: As gas cools due to Brems. radiation, volume of the cone contracts
P_brems = V_2 * P_brems_perV(T_2, rho_N_ArII(T_2, p_ini, n_tot_ini), rho_N_ArII(T_2, p_ini, n_tot_ini), 1)
print("P_brems =", "{:e}".format(P_brems), "W")

    # Aside #1: Graph Brems. radiation with temperature
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

    # Aside #2: Calculate plasma frequency to compare to Brems. light emission
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

#%% Final version of calculations: Pressure with time. At t=0 laser turns on, at t = 10 ms laser turns off. Brems continues after laser off.
start_time = time.time()

# Initialize variables
timestep = 1e-6                # timestep [s]
time_array = np.array([0])      # time [s]
T_ini = 300.0                   # Initial temperature of the chamber (K)
p_ini = 19.91e5                 # Pressure of the chamber (20 bar to Pa)
p = np.array([p_ini])           # pressure [Pa]
T = np.array([T_ini])           # Temperature [K]
T_2_guess = T_ini               # First temperature guess is the initial temp
p_2 = p_ini                     # For step 6.2: Plasma cone is at the same pressure as the surrounding gas
V = np.array([0])               # Volume of plasma cone [m^3]
E_plasma_array = np.array([0])  # Array of energy in the plasma [J]
E_plasma = 0                    # Energy in the plasma [J]
i = 0                           # iteration index
time_end = 10e-3                # end of the laser pulse at 10 ms
P_loss = 0                      # Brems dissipation power [W]
P_laser = 1540                  # laser power [W]
fast_forward_tripped = 0        # Trip so fast forward only happens once on pressure rise
loss = "brems"                  # Switch that takes either "brems" or "blackbody"


while time_array[i] < 0.02:
    if time_array[i] < time_end: # While laser is on
        E_plasma = E_plasma + P_laser * timestep - P_loss * timestep          # Energy in the plasma ("E_plasma") while laser is on
    else: # laser turns off, LSP radiates away energy
        E_plasma = E_plasma - P_loss * timestep          # Energy in the plasma ("E_plasma") goes down as energy is radiated
    
    

    # STEP 6.2: Have XX J of energy to m_plasma, while keeping constant pressure
    T_2 = solve_for_T2(E_plasma, m_plasma, T[i], p[i], T_2_guess, p_2)[0]  # ----TEST IN PROGRESS # Temperature of plasma after energy addition (K)
    T_2_guess = T_2          # Take the last calculated temperature as the guess for the next polynomials

    # STEP 6.3: Volume of the cone contracts; find new volume (V_2) of cone
    n_tot_ini = m_plasma/MW_Ar                                          # Calculate initial number of moles 
    (Gibbs, x, y_Ar, y_ArII, y_e) = MinimizeGibbs_Ar(T_2, p_2)          # Get dissociation x and mole fractions at temp and press of step 2
    n_tot = n_tot_ini*(x + 1)
    V_2 = V_final(p_2, n_tot, T_2)
    
    # STEP 6.4: Calculate pressure increase in chamber due to expansion of gas in plasma core
    p_4 = p_ini * ((V_chamberV1 - V_plasma)/(V_chamberV1 - V_2))**k
    p_2 = p_4           # Plasma cone is at the same pressure as the surrounding gas, for next iteration

    match loss:
        case "brems":
            # STEP 6.5: As gas heats due to Brems. radiation, volume of the cone increases
            P_loss = V_2 * P_brems_perV(T_2, rho_N_ArII(T_2, p_ini, n_tot_ini), rho_N_ArII(T_2, p_ini, n_tot_ini), 1)
        case "blackbody":
            # calculate blackbody radiation loss, as a sanity check upper bound on power loss
            A = area_cone(V_2, r)
            P_loss = constants.Stefan_Boltzmann * T_2**4 *A # blackbody radiation (emissivity e = 1)

    # Printing values for debugging
    print("E_plasma : ", E_plasma)
    print("P_loss : ", P_loss)    
    print("Pressure : ", p_4)
    print("Temperature : ", T_2)

    # Save all values to vectors and increment time + iterator
    time_array = np.append(time_array, time_array[i] + timestep)
    p = np.append(p, p_4)
    T = np.append(T, T_2)               # Store temperature in vector for plotting
    V = np.append(V, V_2)
    E_plasma_array = np.append(E_plasma_array, E_plasma)
    i += 1

    # To speed up calculations, fast forward to laser off when pressure achieves equilibrium
    if np.abs(p[i] - p[i-1]) < 1 and fast_forward_tripped == 0:
        time_array[i] = time_end
        print("Fast-forward to laser off")
        fast_forward_tripped = 1

# Manual zeroing of curves 
time_array += 0.205
p += 500

end_time = time.time()

# Save curves from model for each run
np.save('LSP179_SPRK50_brems_time_2', time_array)
np.save('LSP179_SPRK50_brems_p_2', p)
np.save('LSP179_SPRK50_brems_T_2', T)
np.save('LSP179_SPRK50_brems_V_2', V)
np.save('LSP179_SPRK50_brems_E_plasma_array_2', E_plasma_array)
print(i)
print(end_time-start_time)


#%% Plot figure



# LW41871 = {
#         'psi_conversion': 103.0e-3,  # V/psi
#         'kPa_conversion': 14.94e-3,  # V/kPa 
#     }
#     LW41374 = {
#         'psi_conversion': 100.0e-3,  # V/psi
#         'kPa_conversion': 14.51e-3,  # V/kPa
#     }



# Plot experimental data #1
# plt.figure()
# df = pd.read_csv('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP pressure data 2024-04-05/WaveData19012.csv', header =2, names=['time','Volts']) # Experimental data from LSP178_SPRK49 (V1 100% power pulse shot, 19.91 bar)
# experiment_time = df[["time"]].to_numpy()
# experiment_pressure = df[["Volts"]].to_numpy()* 1e3 / 14.94e-3 + p_ini # Pa, with P_ini added so this is pressure rise
# plt.plot(experiment_time, experiment_pressure, '.')
# plt.plot(time_array, p)
# plt.xlabel('Time [s]')
# plt.ylabel('Pressure [Pa]')
# #plt.xlim([0,20000])
# #plt.ylim([0,5e11])
# plt.title("Pressure change over time during laser pulse, with radiation loss")

# Plot experimental data #2
#LSP179_SPRK50
plt.figure()
df = pd.read_csv('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP pressure data 2024-04-05/WaveData19014.csv', header =2, names=['time','Volts']) # Experimental data from LSP179
experiment_time = df[["time"]].to_numpy()
experiment_pressure = df[["Volts"]].to_numpy()* 1e3 / 14.94e-3 + p_ini # Pa, with P_ini added so this is pressure rise
plt.plot(experiment_time, experiment_pressure, '.')
plt.plot(time_array, p)
plt.xlabel('Time [s]')
plt.ylabel('Pressure [Pa]')
#plt.xlim([0,20000])
#plt.ylim([0,5e11])
plt.title("Pressure change over time during laser pulse, with radiation loss")


plt.legend(['Experimental', '0D model']) # add Brems too


# Plot experimental data #3
#LSP183_SPRK54
# plt.figure()
# df = pd.read_csv('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP pressure data 2024-04-05/WaveData1913.csv', header =2, names=['time','Volts']) # Experimental data from LSP183_SPRK54
# experiment_time = df[["time"]].to_numpy()
# experiment_pressure = df[["Volts"]].to_numpy()* 1e3 / 14.94e-3 + p_ini # Pa, with P_ini added so this is pressure rise
# plt.plot(experiment_time, experiment_pressure, '.')
# plt.plot(time_array, p)
# plt.xlabel('Time [s]')
# plt.ylabel('Pressure [Pa]')
# #plt.xlim([0,20000])
# #plt.ylim([0,5e11])
# plt.title("Pressure change over time during laser pulse, with radiation loss")


# plt.legend(['Experimental', '0D model']) # add Brems too
# %%
