# -*- coding: utf-8 -*-
#!/usr/bin/env python

import os
import numpy as np
from scipy.optimize import leastsq
from scipy import constants 
from scipy import integrate
from scipy.interpolate import splev, splrep

import read_files

# Murnaghan equation of state
def eos_murnaghan(params, vol):
    #From Phys. Rev. B 28, 5480 (1983)#
    E0, B0, Bp, V0 = params 
    E = E0 + B0/Bp * vol * ((V0/vol)**Bp/(Bp-1.0)+1.0) - V0*B0/(Bp-1.0)
    return E

def eos_murnaghan_pressure(params, vol):
    E0, B0, Bp, V0 = params 
    P = (B0/Bp) * ((vol/V0)**(-Bp) - 1)
    return P

def print_params(label, params, fact):
    E0, B0, Bp, V0 = params
    print(label, ": E0 = %f eV" % (E0))
    print(label, ": B0 = %f GPa" % (B0*160.21765)) # eV.A**-3 -> GPa
    print(label, ": Bp = %f" % (Bp))
    print(label, ": V0 = %f angstrom^3" % (V0))
    print(label, ": a0 = %f angstrom" % (V0/fact)**(1.0/3.0))
    print()

#Учет температурной зависимости обменников #
def Murnaghan(fact, sum_conc):

    Data = read_files.read_a_energy(sum_conc)

    # Начальные приближения для V0, E0, B0, Bp из параболы
    a, b, c = np.polyfit(fact*Data[:,0]**3, Data[:,1], 2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    Bp = 4.0

    x0 = [E0, B0, Bp, V0]

    target = lambda params, y, x: y - eos_murnaghan(params, x)
    murn, ier = leastsq(target, x0, args=(Data[:,1], fact*Data[:,0]**3))
    print_params("Murnaghan", murn, fact)

    n = 0
    while(Data[n,0] < (murn[3]/fact)**(1.0/3.0)):
        n = n + 1

    return murn, n


def J_T_accounting(T_start, T_end, T_step, fact, murn, sum_conc, J, d_a):

    Data = read_files.read_a_energy(sum_conc)
    print(Data)
    V_aprox = np.linspace(min(fact*Data[:,0]**3), max(fact*Data[:, 0]**3), 100)
    P_aprox = eos_murnaghan_pressure(murn, V_aprox)*160.21765
    E_aprox = eos_murnaghan(murn, V_aprox)

    Theta0 = input('Enter the Debye temperature:\n')
    Theta0 = float(Theta0)

    T = np.arange(T_start, T_end, T_step) #Температура

    #############  Рассчет параметра Грюнайзера  ##################

    d1 = np.diff(P_aprox, n=1)/np.diff(V_aprox, n=1)
    d2 = np.diff(d1, n=1)/np.diff(V_aprox[1:], n=1)

    gamma = (-2.0/3.0) - (d2 * V_aprox[2:])/(2*d1[1:])
    gamma = sum(gamma) / len(gamma)

    #############  Зависмость температуры Дебая от температуры  ##################

    Theta = ((murn[3]/V_aprox)**gamma) * Theta0

    #############  Рассчет энергии Гельмгольца для каждой температуры  ##################

    # Функция Дебая
    def f(x):
        return x**3.0/(np.exp(x) - 1.0)

    def D(t, Theta_current):
        I = integrate.quad(f, 0.0, float(Theta_current/t))
        return 3.0 * (t/Theta_current)**3.0 * I[0]

    min_T = np.zeros((T.size, 2))
    for i in range(T.size):
        F = np.zeros(V_aprox.size)
        for j in range(V_aprox.size):

            # Вычисление свободной энергии решетки
            F_vib = (constants.R * T[i] * ((9.0/8.0) * (Theta[j]/T[i]) + 3.0*np.log(1.0 - np.exp(-Theta[j]/T[i])) - D(T[i], Theta[j])))/1000.0 #[кДж/моль]

            # Рассчет полной свободной энергии
            F[j] = (E_aprox[j] * 96.32) + F_vib
 
        # Поиск минимума для свободной энергии при данной температуре
        p_T = np.polyfit(V_aprox, F, 2)
        V0_T = -p_T[1] / (2*p_T[0])
        E0_T = p_T[0]*V0_T**2 + p_T[1]*V0_T + p_T[2]

        min_T[i,0] = V0_T
        min_T[i,1] = E0_T

    T_a = []
    spl = splrep((min_T[:,0]/fact)**(1.0/3.0), T)
    for i in range(Data[:,0].size): 
        if Data[i,0] >= (murn[3]/fact)**(1.0/3.0):
            if splev(Data[i,0], spl) > 0.0:
                T_a += [splev(Data[i,0], spl)]
                J_add = read_files.read_J(d_a, True, sorted(os.listdir(os.getcwd()))[i] + '/*JXC.out', True)
                J = np.column_stack((J, J_add))
    T_a = np.array(T_a)

    return J #J_add, T_a
                #J_all = np.column_stack((J_all, J_all_add))
                #J_all_add = []

    

if J_T_accounting == True:
    # Аппроксимация обменников полиномом второй степени (от температуры)#
    J_p_T = np.zeros((J.shape[0], col_add))
    for i in range(J.shape[0]):
        J_p_T[i,:] = np.polyfit(T_a, J[i, 6:6+col_add], 2)

    # Запись коэффициентов полинома вмесо обменников #
    J = np.column_stack((J[:, 0:6], J_p_T))


