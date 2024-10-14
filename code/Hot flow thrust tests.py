# %% LSP 314
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

atm = 1.01325

df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP V2 WinDAQ data/LSP314_V2_Flow3.xlsx', usecols='A, J, K')

time = df['Relative Time']# s
pressure = df['Bar Omega'] + atm # bar
thrust = df['Newton (not good calibration)'] # newton

fig1, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Absolute pressure (bar)', color=color)
ax1.scatter(time, pressure, color=color, marker='.')
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xlim(0,100)
ax1.set_ylim(0,22.5)
ax1.axvline(x = 73, color = 'r', label= 'QCW laser pulse') # line to show laser shot time

ax2 = ax1.twinx()

color = 'tab:blue'

ax2.set_ylabel('Thrust (N)', color=color)  # we already handled the x-label with ax1
ax2.scatter(time, thrust, color=color, marker='.')
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0,0.5)
fig1.legend()

fig1.savefig('../Thesis/assets/4 experiments/LSP314.pdf')

# %% LSP 315

df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP V2 WinDAQ data/LSP315_V2_Flow4.xlsx', usecols='A, J, K')

time = df['Relative Time']# s
pressure = df['Bar Omega'] + atm # bar
thrust = df['Newton (not good calibration)'] # newton

fig2, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Absolute pressure (bar)', color=color)
ax1.scatter(time, pressure, color=color, marker='.')
ax1.tick_params(axis='y', labelcolor=color)
ax1.set_xlim(0,60)
ax1.set_ylim(0,22.5)
ax1.axvline(x = 50, color = 'r', label= 'QCW laser pulse') # line to show laser shot time

ax2 = ax1.twinx()

color = 'tab:blue'

ax2.set_ylabel('Thrust (N)', color=color)  # we already handled the x-label with ax1
ax2.scatter(time, thrust, color=color, marker='.')
ax2.tick_params(axis='y', labelcolor=color)
ax2.set_ylim(0,0.5)
fig2.legend()

fig2.savefig('../Thesis/assets/4 experiments/LSP315.pdf')
# %%
