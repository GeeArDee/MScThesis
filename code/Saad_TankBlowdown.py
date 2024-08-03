# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import data from excel and plot experimental data

start_time = 54.85
stop_time = 57.85

P_atm = 101325e-5   # Pa
P_ini = 20e5 #18.96e5 + P_atm # Pa, Initial pressure in absolute
P_stop = P_atm*1e5/0.528  # Pa, When 

df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP V2 WinDAQ data/LSP315_V2_Flow4.xlsx', usecols='A, K, L')
plt.figure(1)
df_cut = df.loc[(df['Relative Time'] > start_time) & (df['Relative Time'] < stop_time)]
x_exp = df_cut['Relative Time']-start_time # s
y_exp = df_cut['Bar PCB'] # bar
y_exp_Omega = df_cut['Bar Omega']
# plt.errorbar(x_exp, y_exp+P_ini/1e5, yerr=-y_exp*0.01, label='Experimental Data from PCB transducer', fmt='.') # plot in bar
plt.scatter(x_exp-0.02, y_exp_Omega + P_atm , color='black', marker='.', label='Experimental data')
plt.xlabel('Time (s)')
plt.ylabel('Absolute Pressure (bar)')
#plt.title('Saad Tank Blowdown Curve vs Experimental')
#plt.grid(off)

# Initialize arrays and values
P = np.linspace(P_ini, P_stop, num=1000)
V = 9.68*10**-6 # Thruster volume in m^3, ### should includes volume of tubing??
diam_list = [0.20,0.21,0.22]  # throat diameter, in mm
D = np.array(diam_list)*10**-3 # throat diameter, in m
A = np.pi/4*D**2
T = 300 # K
time = np.zeros((len(P),len(A)))
gamma = 1.67
R = 208.13  # J/kg*K

def sol_isentropic(x, V, A, T, gamma, R):
    numerator = -2*V*(x**((1-gamma)/(2*gamma))-1)
    denominator = ((1-gamma)*R*(T)**(1/2)*A*(gamma/R*(2/(gamma+1))**((gamma+1)/(gamma-1)))**(1/2))
    return numerator/denominator

# Isentropic solution first
for i, p in enumerate(P):
    x = p/P_ini
    time[i] = sol_isentropic(x, V, A, T, gamma, R)

plt.figure(1)
plt.plot(time+0.01, P/1e5, label=diam_list) # plot in bar
plt.legend()
#plt.ylim((0,22.5))
#plt.xlim((0,3))

plt.savefig('Thesis/assets/4 experiments/Saad blowdown fit.pdf')

# %%
