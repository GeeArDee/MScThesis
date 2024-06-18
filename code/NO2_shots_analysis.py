# %%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Argon_Control = ['2211','2215','2217']
NO2_0_55_LSP = ['22120','22122','22124']
NO2_0_24_LSP = ['22126','22128']
NO2_0_12_LSP =['22130','22132']




plt.figure()
plt.xlabel('Time (s)')
plt.ylabel('Pressure (Bar)')
plt.title('Comparing the NO2 LSP shots')
plt.grid(which='both')

def get_data_from_index(oscilloscope_index):
    path = 'C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP pressure data 2024-06-08/WaveData' + oscilloscope_index + '.csv'
    df = pd.read_csv(path, skiprows=3,names=['Time', 'Data'])
    # Converting in Bar
    x = df['Time']
    y = df['Data']/1.494
    return x,y

for oscilloscope_index in Argon_Control:
    x, y = get_data_from_index(oscilloscope_index)
    plt.plot(x,y, label='Argon Control', color='xkcd:purple')

for oscilloscope_index in NO2_0_55_LSP:
    x, y = get_data_from_index(oscilloscope_index)
    plt.plot(x,y, label='0.55 Bar NO2', color='xkcd:red')

for oscilloscope_index in NO2_0_24_LSP:
    x, y = get_data_from_index(oscilloscope_index)
    plt.plot(x,y, label='0.24 Bar NO2', color='xkcd:blue')

for oscilloscope_index in NO2_0_12_LSP:
    x, y = get_data_from_index(oscilloscope_index)
    plt.plot(x,y, label='0.12 Bar NO2', color='xkcd:green')
    
plt.legend()

plt.savefig('Thesis/assets/5 results/NO2_shots_analysis.pdf')