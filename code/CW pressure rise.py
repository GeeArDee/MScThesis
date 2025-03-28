# %% LSP 385 - CW!
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

df = pd.read_excel('C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/CW LSP pressure rise WaveData2261.xlsx', usecols='A, B')

time = df['time']# s
pressure = df['uV']/ 14.94e-3 -1.5 +2000# kPa

start = 3.2/1000
end = 88.3/1000

color = 'tab:red'
plt.figure()
plt.xlabel('Time (s)')
plt.ylabel('Pressure (kPa)')
plt.xlim(0,0.4)
plt.ylim(2000,2018)
plt.axvspan(0.20+start, 0.20+end, color='tab:red', alpha=0.5, lw=0, label='LSP lifetime')
plt.legend()
plt.plot(time,pressure, '.', color=color)

plt.savefig('../Thesis/assets/4 experiments/CW pressure rise.pdf')
# %%
