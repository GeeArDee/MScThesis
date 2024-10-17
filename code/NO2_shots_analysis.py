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
plt.ylabel('Dynamic Pressure (bar)')
# plt.title('Comparing the NO2 LSP shots')
plt.grid(which='both')
plt.xlim(0,0.4)
plt.ylim(0,0.16)

def get_data_from_index(oscilloscope_index):
    path = 'C:/Users/gdub5/OneDrive/McGill/Thesis/Experimental things/LSP pressure data 2024-06-08/WaveData' + oscilloscope_index + '.csv'
    df = pd.read_csv(path, skiprows=3,names=['Time', 'Data'])
    # Converting in Bar
    x = df['Time']
    y = df['Data']/1.494
    return x,y



x_NO2_0_55_1, y_NO2_0_55_1 = get_data_from_index('22120')
x_NO2_0_55_2, y_NO2_0_55_2 = get_data_from_index('22122')
x_NO2_0_55_3, y_NO2_0_55_3 = get_data_from_index('22124')

df_concat = pd.concat((y_NO2_0_55_1, y_NO2_0_55_2, y_NO2_0_55_3), axis=1)
y_NO2_0_55__avg = df_concat.mean(axis=1)

plt.plot(x_NO2_0_55_1,y_NO2_0_55__avg, label='0.55 bar $\mathregular{NO_2}$', color='xkcd:dark orange')


x_NO2_0_24_1, y_NO2_0_24_1 = get_data_from_index('22126')
x_NO2_0_24_2, y_NO2_0_24_2 = get_data_from_index('22128')

df_concat = pd.concat((y_NO2_0_24_1, y_NO2_0_24_2), axis=1)
y_NO2_0_24__avg = df_concat.mean(axis=1)

plt.plot(x_NO2_0_24_1,y_NO2_0_24__avg, label='0.24 bar $\mathregular{NO_2}$', color='xkcd:orange')


x_NO2_0_12_1, y_NO2_0_12_1 = get_data_from_index('22130')
x_NO2_0_12_2, y_NO2_0_12_2 = get_data_from_index('22132')

df_concat = pd.concat((y_NO2_0_12_1, y_NO2_0_12_2), axis=1)
y_NO2_0_12__avg = df_concat.mean(axis=1)

plt.plot(x_NO2_0_12_1,y_NO2_0_12__avg, label='0.12 bar $\mathregular{NO_2}$', color='xkcd:light orange')

x_Ar1, y_Ar1 = get_data_from_index('2211')
x_Ar2, y_Ar2 = get_data_from_index('2215')
x_Ar3, y_Ar3 = get_data_from_index('2217')

df_concat = pd.concat((y_Ar1, y_Ar2, y_Ar3), axis=1)
y_Ar_avg = df_concat.mean(axis=1)
plt.plot(x_Ar1,y_Ar_avg, label='Argon control', color='xkcd:purple')


# for oscilloscope_index in NO2_0_55_LSP:
#     x, y = get_data_from_index(oscilloscope_index)
#     plt.plot(x,y, label='0.55 Bar NO2', color='xkcd:red')

# for oscilloscope_index in NO2_0_24_LSP:
#     x, y = get_data_from_index(oscilloscope_index)
#     plt.plot(x,y, label='0.24 Bar NO2', color='xkcd:blue')

# for oscilloscope_index in NO2_0_12_LSP:
#     x, y = get_data_from_index(oscilloscope_index)
#     plt.plot(x,y, label='0.12 Bar NO2', color='xkcd:green')
    
plt.legend()


plt.savefig('../Thesis/assets/4 experiments/NO2_shots_analysis.pdf')
# %%
