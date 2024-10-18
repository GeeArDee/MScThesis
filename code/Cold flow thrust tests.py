# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

atm = 1.01325

df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/Thrust tests/coldflow 5 20 bar -good.xlsx', usecols='A, H, I')

time = df['Relative Time']# s
pressure = df['Bar zeroed'] + atm # bar
thrust = df['Newton zeroed'] # newton

fig1, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Absolute pressure (bar)', color=color)
ax1.scatter(time, pressure, color=color, marker='.')
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xlim(0,20)
ax1.set_ylim(0,22.5)

ax2 = ax1.twinx()

color = 'tab:blue'

ax2.set_ylabel('Thrust (N)', color=color)  # we already handled the x-label with ax1
ax2.scatter(time, thrust, color=color, marker='.')
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0,1.4)

fig1.savefig('../Thesis/assets/4 experiments/Example thrust 20 bar.pdf')

# Making the 35-point thrust vs pressure graph and curve fitting

    # Starting with 5-point one I have:
df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/Thrust tests/_Pressure-Thrust graph.xlsx', usecols='C, D')

pressure_35pts = df['Bar'] + atm # bar
thrust_35pts = df['Newton'] # newton

fig2, ax1 = plt.subplots()
ax1.set_xlabel('Absolute pressure (bar)')
ax1.set_ylabel('Thrust (N)')
ax1.scatter(pressure_35pts, thrust_35pts, marker='.', color='black', label='Experimental data')
ax1.set_xlim(0,40)
ax1.set_ylim(0,2)

# curve fit linear
def func(x, a, b):
    return a*x + b

params = curve_fit(func, pressure_35pts, thrust_35pts)
[a, b] = params[0]
print(a)
print(b)

x = np.linspace(0,40)

ax1.plot(x, func(x,a,b), color='green', label='Linear curve fit', linestyle='dashed')
ax1.legend()

fig2.savefig('../Thesis/assets/4 experiments/pressure-thrust graph.pdf')

# %%
df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/Thrust tests/coldflow 13, 20 bar 5N cell, 200g no cup.xlsx', usecols='A, D')

fig3.savefig('../Thesis/assets/4 experiments/hysterisis graph.pdf')