# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 15:09:47 2023

@author: gdub5
"""

from tabulate import tabulate


def p(m, R_specific, T, V):
    #  p V = m R_s T
    return m*R_specific*T/V


def T(Q, m, c, T_0):
    #  Q = m c delta_T
    return Q/(m*c)+T_0


V_0 = 0.0004  # m^3
MW_Ar = 39.95  # kg/kmol
p_0 = 20e5  # Pa
T_0 = 300  # K
R_u = 8314.462  # J/K-kmol
R_Ar = R_u/MW_Ar
cv = 3/2 * R_Ar
cp = 5/2 * R_Ar
k = cp/cv
Q = 0.8*30  # J

V_A0_list = [1e-6]  # V_0*1/100, V_0*1/10000
V_A0_list_Str = ['1e-6']  # 'V_0*1/100', 'V_0*1/10000'

for i, V_A0 in enumerate(V_A0_list):
    V_B0 = V_0-V_A0
    m_A = p_0*V_A0/(R_Ar*T_0)
    m_tot = p_0*V_0/(R_Ar*T_0)
    m_B = m_tot-m_A

    #  STATE 1
    V_A1 = V_A0
    V_B1 = V_B0
    p_B1 = p_0
    T_A1 = T(Q, m_A, cv, T_0)
    p_A1 = p(m_A, R_Ar, T_A1, V_A0)

    #  STATE 2
    V_A2 = (((p_B1*V_B1**k)/(p_A1*V_A1**k))**(1/k)+1)**(-1)*V_0
    p_A2 = (p_A1*V_A1**k)/(V_A2**k)
    p_B2 = p_A2
    V_B2 = ((p_B1*V_B1**k)/p_B2)**(1/k)
    T_A2 = p_A2*V_A2/(m_A*R_Ar)
    T_B2 = p_B2*V_B2/(m_B*R_Ar)

    # STATE 3
    T_A3 = T(Q, m_A, cp, T_0)
    p_A3 = p_0
    V_A3 = m_A*R_Ar*T_A3/p_A3
    V_B3 = V_0 - V_A3
    p_B3 = p_0*V_B0**k/(V_B3**k)
    T_B3 = p_B3*V_B3/(m_B*R_Ar)

    #  STATE 4
    T_4 = T(Q, m_tot, cv, T_0)
    p_4 = p(m_tot, R_Ar, T_4, V_0)

    #  Printing
    table = [['STATE', 'p [kPa]', 'T [K]'],
             ['A1', p_A1/1000, T_A1],
             ['A2', p_A2/1000, T_A2],
             ['B2', p_B2/1000, T_B2],
             ['A3', p_A3/1000, T_A3],
             ['B3', p_B3/1000, T_B3],
             ['4', p_4/1000, T_4]]
    print(V_A0_list_Str[i])
    print(tabulate(table, headers='firstrow', tablefmt='fancy_grid'))
