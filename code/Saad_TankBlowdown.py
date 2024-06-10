# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import data from excel and plot experimental data

start_time = 54.85
stop_time = 55.38

df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP V2 WinDAQ data/LSP315_V2_Flow4.xlsx', usecols='A, L')
plt.figure(1)
df_cut = df.loc[(df['Relative Time'] > start_time) & (df['Relative Time'] < stop_time)]
x_exp = df_cut['Relative Time']-start_time # s
y_exp = df_cut['Bar PCB'] # bar
plt.errorbar(x_exp, y_exp, yerr=-y_exp*0.01, label='Experimental Data from PCB transducer', fmt='.') # plot in bar cool
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Bar)')
plt.title('Saad Tank Blowdown Curve vs Experimental')
plt.grid(which='both')

# Initialize arrays and values
start_press = 0  # Bar
stop_press = -7 # Bar
P = np.linspace(start_press, stop_press, num=1000)*1e5
P_atm = 101325   # Pa
P_ini = 18.96*1e5 + P_atm # Pa, Initial pressure in absolute
V = 9.68*10**-6 # Thruster volume in m^3, ### should includes volume of tubing??

diam_list = [0.2,0.3,0.35]  # throat diameter, in mm

D = np.array(diam_list)*10**-3 # throat diameter, in m
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
    time[i] = sol_isentropic(x, V, A, T, gamma, R)#+start_time+0.03

plt.figure(1)
plt.plot(time, P/1e5 - 0.52664, label=diam_list) # plot in bar
plt.legend()

# %%
