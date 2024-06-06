# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import data from excel and plot experimental data
start_time = 54.89
stop_time = 55.36
df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP V2 WinDAQ data/LSP315_V2_Flow4.xlsx', usecols='A, L')
plt.figure(1)
df_cut = df.loc[(df['Relative Time'] > start_time) & (df['Relative Time'] < stop_time)]
x_exp = df_cut['Relative Time']
y_exp = df_cut['Bar PCB']*1e5
plt.plot(x_exp, y_exp)

# Initialize arrays and values
start_press = 0  # Bar
stop_press = -14 # Bar
P = np.linspace(start_press, stop_press, num=1000)*1e5
P_atm = 101325   # Pa
P_ini = 18.96*1e5  # Pa
V = 9.68*10**-6 # Thruster volume in m^3, ### should includes volume of tubing??
D = np.array([0.7, 1])*10**-3 # throat diameter, in m
A = np.pi*(D/2)**2
T = 300 # K
time = np.zeros((len(P),len(A)))

gamma = 1.67
R = 208.13  # J/kg*K

def sol_isentropic(x, V, A, T, gamma, R):
    return -2*V*(x**((1-gamma)/2*gamma)-1)/((1-gamma)*R*np.sqrt(T)*A*np.sqrt(gamma/R*(2/(gamma+1))**((gamma+1)/(gamma-1))))

# Isentropic solution first

for i, p in enumerate(P):
    P_res = P_ini + p
    x = P_res/P_ini
    time[i] = sol_isentropic(x, V, A, T, gamma, R)+start_time

plt.figure(1)
plt.plot(time, P)


# %%
