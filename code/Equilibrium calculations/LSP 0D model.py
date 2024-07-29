#%%
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 14:29:59 2023

@author: gdub5
"""

import matplotlib.pyplot as plt
import numpy as np
from CEA_LSP_polynomials import CEA_LSP_polynomials
from scipy.optimize import minimize, Bounds
import pandas as pd
import sympy as sp
from labellines import labelLines

# Physical constants

n_pts = 10000  # Number of points for curves
R = 8314.462  # Universal gas constant [J/kmol-K]
T_ref = 298.15  # K
p_ref = 100000  # 1 bar = 1e5 Pa
MW_Ar = 39.948
MW_ArII = 39.9474514
MW_e = 0.000548579903

# Extract properties from polynomials


def cp(compound, T):
    return CEA_LSP_polynomials(compound, T)[0]*R


def h_bar_i(compound, T):
    return CEA_LSP_polynomials(compound, T)[1]*R*T


def s_bar_0_i(compound, T):
    return CEA_LSP_polynomials(compound, T)[2]*R


def entropy_bar_i(compound, T, p_i):
    return s_bar_0_i(compound, T) - R*np.log(p_i/p_ref, where=p_i > 0)


def gibbs_bar_i(compound, T, p_i):
    return h_bar_i(compound, T) - T*entropy_bar_i(compound, T, p_i)


def gibbs_bar_mixture_H2(x, p, T):
    n_H = 2*x
    n_H2 = 1-x
    n_tot = n_H + n_H2

    y_H = n_H/n_tot
    y_H2 = n_H2/n_tot

    p_H = y_H*p
    p_H2 = y_H2*p

    return (n_H*gibbs_bar_i("H", T, p_H) + n_H2*gibbs_bar_i("H2", T, p_H2))/1e6


def gibbs_bar_mixture_Ar(x, p, T):
    n_Ar = 1-x
    n_ArII = x
    n_e = x
    n_tot = n_Ar + n_ArII + n_e

    y_Ar = n_Ar/n_tot
    y_ArII = n_ArII/n_tot
    y_e = n_e/n_tot

    p_Ar = y_Ar*p
    p_ArII = y_ArII*p
    p_e = y_e*p

    gibbs_mix = (n_Ar*gibbs_bar_i("Ar", T, p_Ar) + n_ArII *
                 gibbs_bar_i("Arp", T, p_ArII) + n_e*gibbs_bar_i('e', T, p_e))/1e6
    return gibbs_mix

# def K(T, p):
#     return

# def alpha(T, p):
#     return np.sqrt(K(T))/np.sqrt(K(t) + p/pref)

# # 'compound' values: 'H2', 'H', 'Hp', 'e', 'Ar', 'Arp'
# test = gibbs_bar_i('Ar', 4000, 10e5)
# print(test/1e6)  # Gibbs in MJ

# %% 1) Check that NASA polynomials work to give c_p(T), h(T), and s0(T) for H2


T = np.arange(500, 20001, 500, dtype=float)  # temperature in K

# compound = ['H2', 'H']  # , 'H', 'H2', 'Hp', 'e', 'Ar', 'Arp'

size = np.size(T)
cp_values = np.zeros(size)
h_values = np.zeros(size)
s0_values = np.zeros(size)

compound = 'H2'
T = np.arange(500, 20001, 500, dtype=float)  # temperature in K

for i, t in enumerate(T):
    cp_values[i] = cp(compound, t)
    h_values[i] = h_bar_i(compound, t)
    s0_values[i] = s_bar_0_i(compound, t)

data = np.array([cp_values, h_values/1000, s0_values]).T
df = pd.DataFrame(data, columns=['c_p', 'h', 's0'], index=T)

# display dataframe
print(compound)
print(df)

# %% 2) Table of Gibbs function (per kmol) for diatomic hydrogen

P = np.array([0.1e5, 1e5, 10e5])

s_H2 = np.zeros((np.size(T), np.size(P)))
g_H2 = np.zeros((np.size(T), np.size(P)))

for i, t in enumerate(T):
    for j, p in enumerate(P):
        g_H2[i, j] = gibbs_bar_i(compound, t, p)  # Gibbs function g(T,p)

df_H2 = pd.DataFrame(
    g_H2/1e6, columns=['g_H2, 0.1 bar', '1 bar', '10 bar'], index=T)

# %% 3) Hydrogen in a piston problem
# Draw plot
bounds = Bounds(lb=0.0, ub=1.0)
x = np.linspace(0, 1, 10000)
T = 4000
P = np.array([0.1e5, 1e5, 10e5])  # pressure in pa

for p in P:
    G = gibbs_bar_mixture_H2(x, p, T)
    plt.plot(x, G)  # G in MJ/kmol
    xi = minimize(gibbs_bar_mixture_H2, 0.5, args=(p, T), bounds=bounds).x[0]
    print(xi)
    print(f'Mole Fraction of H at {p/1e5} bar: {2*xi/(1+xi)}')
    print(f'Mole Fraction of H2 at {p/1e5} bar: {(1-xi)/(1+xi)}\n')

plt.xlim(0, 1)
plt.legend(['0.1 bar', '1 bar', '10 bar'])
plt.xlabel('x (degree of dissociation)')
plt.ylabel('G (MJ)')
plt.savefig('../../Thesis/assets/2 models/Gibbs.pdf')

# %% 4) Solve, minimizing gibbs energy with Argon

# x = np.linspace(0, 1, 10000)
T = [15000, 16000]
P = [0.1e5]

# G = gibbs_bar_mixture_Ar(x, P, T)
# plt.plot(x, G)  # G in MJ/kmol


def MinimizeGibbs_Ar(T, P):
    bounds = Bounds(lb=0.0, ub=1.0)
    solution = minimize(gibbs_bar_mixture_Ar, 0.5, args=(P, T), bounds=bounds)
    Gibbs = solution.fun
    x = solution.x[0]
    yAr = (1-x)/(x+1)
    yArII = x/(x+1)
    ye = x/(x+1)
    return Gibbs, yAr, yArII, ye


# # Verify if mole fractions are the same as Higgins'


# g_Ar = np.zeros((np.size(T), np.size(P)))
# yAr = np.zeros((np.size(T), np.size(P)))
# yArII = np.zeros((np.size(T), np.size(P)))
# ye = np.zeros((np.size(T), np.size(P)))

# for i, t in enumerate(T):
#     for j, p in enumerate(P):
#         [g_Ar[i, j], yAr[i, j], yArII[i, j], ye[i, j]
#          ] = MinimizeGibbs_Ar(t, p)  # Gibbs function G(T,P)


# Define these functions to be able to calculate cp = dh/dT (done with sympy)


def MW_avg_Ar(T, P):
    (Gibbs, yAr, yArII, ye) = MinimizeGibbs_Ar(T, P)
    return yAr*MW_Ar + yArII*MW_ArII + ye*MW_e


def h_bar_i_symbolic(compound, T):
    t = sp.symbols('t')
    return CEA_LSP_polynomials(compound, T)[4]*R*t


def h_mix_bar_Ar(T, P):
    (Gibbs, yAr, yArII, ye) = MinimizeGibbs_Ar(T, P)
    return yAr*h_bar_i_symbolic('Ar', T) + yArII*h_bar_i_symbolic('Arp', T) + ye*h_bar_i_symbolic('e', T)


def h_mix_mass_Ar(T, P):
    h_expr = h_mix_bar_Ar(T, P)/MW_avg_Ar(T, P)
    return h_expr


def cp_mix_mass(T, P):
    t = sp.symbols('t')
    expr = sp.diff(h_mix_mass_Ar(T, P), t)
    collected_expr = sp.collect(expr, t)
    return collected_expr/1e3


def cp_mix_mass_take2(T, P):  # Cp in J/(kg*K)
    t = sp.symbols('t')
    for i, temp in enumerate(T):
        for j, press in enumerate(P):
            h = sp.lambdify(t, h_mix_mass_Ar(temp, press), "numpy")
            h_results[i, j] = h(temp)
    return np.gradient(h_results, T, axis=0)


# temperature in K, change 1000 to 100 for smooth curve
T = np.arange(1000, 20001, 1000, dtype=float)
P = np.array([0.1e5, 1e5, 10e5])  # pressure in Pa
h_results = np.zeros((np.size(T), np.size(P)))
cp_results = np.zeros((np.size(T), np.size(P)))

t = sp.symbols('t')

for j, press in enumerate(P):
    for i, temp in enumerate(T):
        h = sp.lambdify(t, h_mix_mass_Ar(temp, press), "numpy")
        h_results[i, j] = h(temp)
        cp = sp.lambdify(t, cp_mix_mass(temp, press), "numpy")
        cp_results[i, j] = cp(temp)

cp_take2 = cp_mix_mass_take2(T, P)
CEA_temp = np.arange(1000, 19001, 1000, dtype=float)
CEA_cp = (np.array([[0.5203,   0.5203,   0.5203,   0.5203,
                    0.5204, 0.5291, 0.5913, 0.8588, 1.6597, 3.4980,
                    6.8383, 10.9168,  11.8562, 8.0516, 4.2650,  2.3774,
                    1.6076,  1.3004,  1.1754],
                    [0.5203, 0.5203,  0.5203, 0.5203, 0.5205, 0.5231, 0.5428,
                    0.6270, 0.8814,   1.4692,   2.5990, 4.4561,
                    6.9595,  9.2153, 9.5039, 7.4502, 4.8822, 3.0987, 2.1205],
                    [0.5203,  0.5203, 0.5203, 0.5203, 0.5204,  0.5212,  0.5274,
                    0.5537,  0.6351,  0.8223,  1.1837, 1.7966, 2.7290,  4.0008,
                    5.5029,  6.8846,  7.5734,  7.1748,  5.9329]]).T)*1e3
Higgins_temp = np.arange(1000, 20001, 1000, dtype=float)
Higgins_cp = (np.array([[0.52033093,
                        0.520330931,
                        0.520330932,
                        0.520336321,
                        0.520789321,
                        0.529063027,
                        0.591254816,
                        0.858763504,
                        1.659684251,
                        3.49799123,
                        6.838266807,
                        10.91674465,
                        11.85609881,
                        8.051556419,
                        4.264955785,
                        2.377436027,
                        1.607562186,
                        1.3003758,
                        1.175405444,
                        1.125748033], [0.52033093,
                                       0.520330931,
                                       0.52033093,
                                       0.520332639,
                                       0.520475881,
                                       0.523092435,
                                       0.542752545,
                                       0.626960739,
                                       0.881431,
                                       1.469147738,
                                       2.599002989,
                                       4.456035764,
                                       6.959441245,
                                       9.215302524,
                                       9.503830062,
                                       7.45018046,
                                       4.882157917,
                                       3.098670678,
                                       2.120477833,
                                       1.627078403], [0.52033093,
                                                      0.520330931,
                                                      0.520330929,
                                                      0.520331475,
                                                      0.520376763,
                                                      0.521204368,
                                                      0.527414761,
                                                      0.553654786,
                                                      0.63512524,
                                                      0.82234088,
                                                      1.183741816,
                                                      1.796625672,
                                                      2.728958959,
                                                      4.000749598,
                                                      5.502839642,
                                                      6.884522833,
                                                      7.573396024,
                                                      7.174766799,
                                                      5.932919521,
                                                      4.496251534
                                                      ]]).T)*1e3

# plt.figure()
# plt.plot(T, h_results)

# plt.figure()
# plt.plot(T, cp_results)

plt.figure()
plt.plot(T, cp_take2)

plt.xlabel('Temperature (K)')
plt.ylabel('Cp')


# Add Higgins and CEA datapoints to plot
plt.plot(CEA_temp, CEA_cp, 'o')

plt.legend(['0.1 bar', '1 bar', '10 bar', 'CEA data'])
# plt.plot(Higgins_temp, Higgins_cp, 'o')

# UNCOMMENT TO OVERWRITE FANCY SMOOTH FIG
#plt.savefig('../../Thesis/assets/2 models/Cp_compare.pdf')

# %% 5)Test problems for heat addition part 2, variable cp


# %% 6) 0D Heat balance model

# Initial values
T_plasma = 15000  # K
P_test = 20e5  # 20 bar to Pa
Lzr_power = 3000  # W
Lzr_duration = 10/1000  # s
delta_q = Lzr_power*Lzr_duration  # J

# Knowing cp from Equilibrium, find mass and volume of plasma ball, assuming constant pressure

# Calculate rate of heat loss via bremsstralung.
# Is plasma ball optically thick or thin? -> volumetric/surface emitter
# Rate of heat loss? Suppose it's thick (a surface emitter)?


# delta_q = m*(h2-h1)

# Volume, then mass of plasma ball
d = 2/100  # m
#V = np.pi*
