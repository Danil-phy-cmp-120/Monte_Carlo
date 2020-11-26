# -*- coding: utf-8 -*-
#!/usr/bin/env python

import random
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import simps
from scipy.optimize import leastsq
from scipy.interpolate import splev, splrep
import math
from scipy import constants 
from scipy import integrate

##########################################################
L = 6
d_a = 2.0

T_start = 0.0
T_end = 600.0
T_step = 20.0

nmcs = 10000
nstep = 400
mcs_start = int(nmcs/10)

J_T_accounting = False
##########################################################

R = 8.314459848  # Дж/(моль*K) 
kB = 0.086173303 # мэВ / К

flag1 = False
flag2 = False
flag3 = False
flag4 = False
Data_lattice = []
J_all = []
J = []
IQ = []
conc = []
atoms = []
atoms_cpa = []


if J_T_accounting == False:
    f = open('*JXC.out',"r")
    lines = f.readlines()
    f.close()

    # Считываем вектора трансляций, базис #
    for line in lines:
        inp = line.split()
        if len(inp) == 0:
            continue
        if flag1 == True and len(inp) == 5 and inp[0] == '(':
            Data_lattice += [float(inp[1].replace(',','')), float(inp[2].replace(',','')), float(inp[3])]
        if inp[0] == '<INIT_MOD_LATTICE>' and len(inp) >= 1:
            flag1 = True
        if inp[len(inp)-1] == '2*pi/a' and len(inp) > 1:
            break

    primitive_vectors = np.array(Data_lattice[:9])
    primitive_vectors.shape = (3, 3)
    #print(primitive_vectors)

    A = np.linalg.inv(primitive_vectors) # Матрица перехода от декартового базиса к текущему
    #print(A)
 
    basis = np.array(Data_lattice[9:len(Data_lattice)])
    basis.shape = (int(len(basis)/3), 3)
    #print(basis)

    # Считываем тип атомов #
    for line in lines:
        inp = line.split()
        if flag2 == True  and len(inp) >= 9:
            for i in range(int(inp[len(inp)-1]) - int(inp[8]) + 1):
                IQ += [inp[1]]
                conc += [float(inp[7])]
        if len(inp) > 3:
            if inp[0] == 'type' and inp[1] == 'TXTT' and inp[2] == 'NL':
                flag2 = True
        if len(inp) == 1:
            if inp[0] == 'mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm':
                break
    conc = np.array(conc)
    print(IQ)


    # Считываем обменные интегралы #
    for line in lines:
        inp = line.split()
        if len(inp) == 0:
            continue
        if flag3 == True and len(inp) > 10:
            J_all += ["{:.3f}".format(float(inp[5])), "{:.3f}".format(float(inp[6])), "{:.3f}".format(float(inp[7])), "{:.3f}".format(float(inp[8])), "{:.8f}".format(float(inp[10]))]
            flag3 = False
        if inp[0] == 'ITAUIJ' and inp[1] == 'ITAUJI' and len(inp) >= 1:
            flag3 = True
        if inp[0] == 'IQ' and inp[6] == 'JQ' and len(inp) == 12:
            J_all += ["{:.0f}".format(float(inp[5])), "{:.0f}".format(float(inp[11]))]

    J_all = np.array(J_all, dtype=float)
    J_all.shape = (int(J_all.size/7), 7)
    #print(J_all)

if J_T_accounting == True:
    #Учет температурной зависимости обменников #
    # Murnaghan equation of state
    def eos_murnaghan(params, vol):
        'From Phys. Rev. B 28, 5480 (1983)'
        E0, B0, Bp, V0 = params 
        E = E0 + B0/Bp * vol * ((V0/vol)**Bp/(Bp-1.0)+1.0) - V0*B0/(Bp-1.0)
        return E

    def eos_murnaghan_pressure(params, vol):
        E0, B0, Bp, V0 = params 
        P = (B0/Bp) * ((vol/V0)**(-Bp) - 1)
        return P

    def print_params(label, params):
        E0, B0, Bp, V0 = params
        print(label, ": E0 = %f eV" % (E0))
        print(label, ": B0 = %f GPa" % (B0*160.21765)) #eV.A**-3 -> GPa
        print(label, ": Bp = %f" % (Bp))
        print(label, ": V0 = %f angstrom^3" % (V0))
        print(label, ": a0 = %f angstrom" % (V0/fact)**(1.0/3.0))
        print()


    #############  Настройки  ##################
    Basis = 'fcc'
    Theta0 = 340.5926 #Температура Дебая
    T = np.arange(T_start, T_end, T_step) #Температура
    ############################################

    #############  Считывание параметра решетки и энергии  ##################
    Data = []
    path = os.listdir(os.getcwd())
    for d in range(len(path)):

        if os.path.isdir(path[d]) == True and path[d] != 'input':

            f = open(path[d] + '/*SCF.out',"r")
            lines = f.readlines()
            f.close()

            for line in lines:
                inp = line.split()
                if len(inp) > 3 and inp[0] == 'ETOT' and inp[2] == 'SCF':
                    E = inp[1]

            Data += [float(path[d]), (float(E)*13.605692)/4]

    Data = np.array(Data)
    Data.shape = ((Data.size/2,2))
    Data = Data[Data[:,0].argsort()]
    #print(Data)


    if Basis == "sc": fact = 1.0
    if Basis == "bcc": fact = 0.5
    if Basis == "fcc": fact = 0.25

    # Начальные приближения для V0, E0, B0, Bp из параболы
    a, b, c = np.polyfit(fact*Data[:,0]**3, Data[:,1], 2)
    V0 = -b/(2*a)
    E0 = a*V0**2 + b*V0 + c
    B0 = 2*a*V0
    Bp = 4.0

    x0 = [E0, B0, Bp, V0]

    target = lambda params, y, x: y - eos_murnaghan(params, x)
    murn, ier = leastsq(target, x0, args=(Data[:,1], fact*Data[:,0]**3))
    print_params("Murnaghan", murn)

    V_aprox = np.linspace(min(fact*Data[:,0]**3), max(fact*Data[:, 0]**3), 100)
    P_aprox = eos_murnaghan_pressure(murn, V_aprox)*160.21765
    E_aprox = eos_murnaghan(murn, V_aprox)
    #print(P_aprox)

    #V_aprox = fact*Data[:,0]**3
    #P_aprox = eos_murnaghan_pressure(murn, V_aprox)*160.21765

    #############  Рассчет параметра Грюнайзера  ##################

    d1 = np.diff(P_aprox, n=1)/np.diff(V_aprox, n=1)
    d2 = np.diff(d1, n=1)/np.diff(V_aprox[1:], n=1)

    gamma = (-2.0/3.0) - (d2 * V_aprox[2:])/(2*d1[1:])
    gamma = sum(gamma) / len(gamma)

    #############  Зависмость температуры Дебая от температуры  ##################

    Theta = ((x0[3]/V_aprox)**gamma) * Theta0

    #############  Рассчет энергии Гельмгольца для каждой температуры  ##################

    # Функция Дебая
    def f(x):
        return x**3.0/(np.exp(x) - 1.0)

    def D(t, Theta_current):
        I = integrate.quad(f, 0.0, float(Theta_current/t))
        return 3.0 * (t/Theta_current)**3.0 * I[0]


    fig = plt.figure()
    ax = fig.add_subplot(111)

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

        ax.plot(V_aprox, F)


    ax.plot(min_T[:,0], min_T[:,1], 'ro')

    ax.grid()
    plt.tight_layout()
    #ax.set_xlim(0.0,500.0)
    #ax.set_ylim(0.0,3.0)
    fig.savefig('F(V).png', dpi=300)


    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot((min_T[:,0]/fact)**(1.0/3.0), T, '-ro')
    ax.grid()
    plt.tight_layout()
    fig.savefig('a(T).png', dpi=300)



    J_all_add = []
    col_add = 0
    T_a = []
    spl = splrep((min_T[:,0]/fact)**(1.0/3.0), T)
    for i in range(Data[:,0].size): 

        if Data[i,0] >= (murn[3]/fact)**(1.0/3.0) and flag4 == True:

            if os.path.isdir(str(Data[i,0])) == True:

                T_a += [splev(Data[i,0], spl)]
                col_add = col_add + 1

                f = open(str(Data[i,0]) + '/*JXC.out',"r")
                lines = f.readlines()
                f.close()

                for line in lines:
                    inp = line.split()
                    if len(inp) == 0:
                        continue
                    if flag3 == True and len(inp) > 10:
                        J_all_add += ["{:.8f}".format(float(inp[10]))]
                        flag3 = False
                    if inp[0] == 'ITAUIJ' and inp[1] == 'ITAUJI' and len(inp) >= 1:
                        flag3 = True

                J_all_add = np.array(J_all_add, dtype=np.float)
                J_all = np.column_stack((J_all, J_all_add))
                J_all_add = []
                



        if Data[i,0] >= (murn[3]/fact)**(1.0/3.0) and flag4 == False:

            if os.path.isdir(str(Data[i,0])) == True:

                flag4 = True
                T_a += [splev(Data[i,0], spl)]
                col_add = col_add + 1

                first_point = Data[i,0]

                f = open(str(Data[i,0]) + '/*JXC.out',"r")
                lines = f.readlines()
                f.close()

                # Считываем вектора трансляций, базис #
                for line in lines:
                    inp = line.split()
                    if len(inp) == 0:
                        continue
                    if flag1 == True and len(inp) == 5 and inp[0] == '(':
                        Data_lattice += [float(inp[1].replace(',','')), float(inp[2].replace(',','')), float(inp[3])]
                    if inp[0] == '<INIT_MOD_LATTICE>' and len(inp) >= 1:
                        flag1 = True
                    if inp[len(inp)-1] == '2*pi/a' and len(inp) > 1:
                        break

                primitive_vectors = np.array(Data_lattice[:9])
                primitive_vectors.shape = (3, 3)
                #print(primitive_vectors)

                A = np.linalg.inv(primitive_vectors) # Матрица перехода от декартового базиса к текущему
                #print(A)
 
                basis = np.array(Data_lattice[9:len(Data_lattice)])
                basis.shape = (len(basis)/3, 3)
                #print(basis)

                # Считываем тип атомов #
                for line in lines:
                    inp = line.split()
                    if flag2 == True  and len(inp) >= 9:
                        for i in range(int(inp[len(inp)-1]) - int(inp[8]) + 1):
                            IQ += [inp[1]]
                            conc += [float(inp[7])]
                    if len(inp) > 3:
                        if inp[0] == 'type' and inp[1] == 'TXTT' and inp[2] == 'NL':
                            flag2 = True
                    if len(inp) == 1:
                        if inp[0] == 'mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm':
                            break
                conc = np.array(conc)
                #print(IQ)


                # Считываем обменные интегралы #
                for line in lines:
                    inp = line.split()
                    if len(inp) == 0:
                        continue
                    if flag3 == True and len(inp) > 10:
                        J_all += ["{:.3f}".format(float(inp[5])), "{:.3f}".format(float(inp[6])), "{:.3f}".format(float(inp[7])), "{:.3f}".format(float(inp[8])), "{:.8f}".format(float(inp[10]))]
                        flag3 = False
                    if inp[0] == 'ITAUIJ' and inp[1] == 'ITAUJI' and len(inp) >= 1:
                        flag3 = True
                    if inp[0] == 'IQ' and inp[6] == 'JQ' and len(inp) == 12:
                        J_all += ["{:.0f}".format(float(inp[5])), "{:.0f}".format(float(inp[11]))]

                J_all = np.array(J_all, dtype=float)
                J_all.shape = (J_all.size/7, 7)
                #print J_all

                


    T_a = np.array(T_a)
   #print(col_add, T_a)

#Учитываем обменники до определенного растояния #
for j in range(J_all[:,5].size): 
     print(J_all[j,5], d_a) 
     if J_all[j,5] <= d_a:
         if J_T_accounting == False:
             J += [float(J_all[j,0]), float(J_all[j,1]), float(J_all[j,2]), float(J_all[j,3]), float(J_all[j,4]), float(J_all[j,5]), float(J_all[j,6])*1000]
         if J_T_accounting == True:
             J += [float(J_all[j,0]), float(J_all[j,1]), float(J_all[j,2]), float(J_all[j,3]), float(J_all[j,4]), float(J_all[j,5])]
             for k in range(col_add):
                 J += [float(J_all[j,6+k])*1000]

J = np.array(J)
if J_T_accounting == True:
    J.shape = (J.size/(6+col_add), (6+col_add))
else:
        J.shape = (int(J.size/7), 7)

#Убираем маленькие обменники #
new_J = []
for i in range(J.shape[0]):
    if abs(J[i,6]) > 0.01:
        if J_T_accounting == False:
            new_J += [J[i,0], J[i,1], J[i,2], J[i,3], J[i,4], J[i,5], J[i,6]]
        if J_T_accounting == True:
            new_J += [J[i,0], J[i,1], J[i,2], J[i,3], J[i,4], J[i,5]]
            for k in range(col_add):
                new_J += [float(J[i,6+k])]
J = np.array(new_J)
if J_T_accounting == True:
    J.shape = (J.size/(6+col_add), (6+col_add))
else:
        J.shape = (int(J.size/7), 7)

J_reverse = J.copy()
for i in range(J_reverse.shape[0]):
    J_reverse[i][0], J_reverse[i][1] = J_reverse[i][1], J_reverse[i][0]
J = np.vstack((J, J_reverse))

J_string = []
new_J = []
for i in range(J.shape[0]):
    j_string = ''
    for j in range(J.shape[1]):
        j_string = j_string + str(J[i,j])

    if j_string not in J_string:
        J_string += [j_string]
        if J_T_accounting == False:
            new_J += [float(J[i,0]), float(J[i,1]), float(J[i,2]), float(J[i,3]), float(J[i,4]), float(J[i,5]), float(J[i,6])]
        if J_T_accounting == True:
            new_J += [float(J[i,0]), float(J[i,1]), float(J[i,2]), float(J[i,3]), float(J[i,4]), float(J[i,5])]
            for k in range(col_add):
                new_J += [float(J[i,6+k])]

J = np.array(new_J)
if J_T_accounting == True:
    J.shape = (J.size/(6+col_add), (6+col_add))
else:
        J.shape = (int(J.size/7), 7)
#print(J)
np.savetxt('J.dat', J, fmt='%.5f')

if J_T_accounting == True:
    # Аппроксимация обменников полиномом второй степени (от температуры)#
    J_p_T = np.zeros((J.shape[0], col_add))
    for i in range(J.shape[0]):
        J_p_T[i,:] = np.polyfit(T_a, J[i, 6:6+col_add], 2)

    # Запись коэффициентов полинома вмесо обменников #
    J = np.column_stack((J[:, 0:6], J_p_T))


for i in range(J.shape[0]): # Перевод векторов из Декатового базиса в базис решетки
    J[i,2:5] = A.dot(J[i,2:5])


for i in range(len(conc*L**3)): # Небольшое изменение концентрации под текущую ячейку чтобы число атомов было int
    conc[i] = round(conc[i]*L**3)/(L**3)
print(conc)


# Транслируем атомы базиса по векторам трансляций #
Lattice = np.zeros((int(sum(conc))*L**3, 7))
atom_index = 1
for j in range(len(IQ)):
    if conc[j] != 1.0:
        for m in range(int(L**3*conc[j])):
            atoms_cpa += [atom_index]
    if conc[j] == 1.0 or j == len(IQ)-1:
        random.shuffle(atoms_cpa) #премешивание атомов на одной позиции
        atoms = atoms + atoms_cpa
        atoms_cpa = []
        for m in range(L**3):
            atoms += [atom_index]

    if j < len(IQ)-1:
        if IQ[j] != IQ[j+1]:
            atom_index = atom_index + 1

#print(atoms)
                    
random.seed() 
n = 0
for i in range(basis[:,0].size):
    for n1 in range(L):
        for n2 in range(L):
            for n3 in range(L):
                C = np.dot(primitive_vectors, np.array([n1, n2, n3]))
                Lattice[n,0] = atoms[n]
                
                Lattice[n,1] = basis[i,0] + C[0]
                Lattice[n,2] = basis[i,1] + C[1]
                Lattice[n,3] = basis[i,2] + C[2]

                # Начальные условия (спин по z)
                #Lattice[n,4] = random.random() 
                #Lattice[n,5] = random.random() 
                #Lattice[n,6] = random.random()
 
                #Lattice[n,4] = 1 
                #Lattice[n,5] = 1 
                Lattice[n,6] = 1 
                n = n + 1


Lattice = sorted(Lattice, key=lambda a_entry: a_entry[0])
Lattice = np.array(Lattice)


for i in range(Lattice.shape[0]): # Перевод векторов из Декатового базиса в базис решетки
    Lattice[i,1:4] = A.dot(Lattice[i,1:4])


# Считывание магнитных моментов из SCF.out #
if J_T_accounting == False:
    f = open('*SCF.out',"r")
    lines = f.readlines()
    f.close()
else:
    f = open(str(first_point) + '/*SCF.out',"r")
    lines = f.readlines()
    f.close()

magmom = np.zeros(len(IQ))
flag = False
for i in range(len(IQ)):
    for line in lines[len(lines)-500:len(lines)]:
        inp = line.split()
        if len(inp) > 5:
            if inp[1] == 'E=' and inp[4] == 'IT=' and inp[6] == IQ[i]:
                flag = True
        if len(inp) > 9:
            if flag == True and inp[0] == 'sum':
                magmom[i] = inp[4]
                flag = False
#print(magmom)


# Сохранение файла задачи для кода C++#
if os.path.exists('input') == False:
    os.mkdir('input')
file = open('input/Settings.dat', "w")
file.write(str(T_start) + ' ')
file.write(str(T_end) + ' ')
file.write(str(T_step) + ' ')
file.write(str(nmcs) + ' ')
file.write(str(nstep) + ' ')
file.write(str(mcs_start) + ' ')

file.write(str(L) + ' ')
file.write(str(int(sum(conc))) + ' ')

file.write(str(J.shape[0]) + ' ') #Число учитываемых обменных интеграллов
file.write(str(J.shape[1])) #Число учитываемых обменных интеграллов

file.close()


file = open('input/Lattice.dat', "wb")
np.savetxt(file, Lattice, fmt='%8f') #Данные о ячейке и спинах атомов
file.close()

file = open('input/J.dat', "wb")
np.savetxt(file, J, fmt='%8f') 
file.close()

file = open('input/J.dat', "wb")
np.savetxt(file, J, fmt='%8f') 
file.close()

file = open('input/magmom.dat', "wb")
np.savetxt(file, magmom, fmt='%8f') 
file.close()


'''
# Сохранение структуры в формате POSCAR #
file = open('POSCAR_Monte_Carlo', "wb")
file.write('Structure_for_Monte-Carlo_calculation' + '\n')
file.write('15.0000000' + '\n')
np.savetxt(file, A.dot(primitive_vectors), fmt='%8f') 
for i in IQ:
   file.write(i + '  ')
file.write('\n')

m = 20
for n in range(len(IQ)):
    for i in Lattice[:,0]:
        if n+1 == int(i):
            m = m + 1
    file.write(str(m) + '  ')
    m = 0
file.write('\n')
file.write('Direct')
file.write('\n')
np.savetxt(file, Lattice[:,1:4]/max(Lattice[:,1]), fmt='%8f') 
file.close()
'''

'''
All_Data = Metropolis.Metropolis_cpp(Lattice, Lattice.shape[0], J, J.shape[0], int(sum(conc)), L, T_start, T_end, T_step, nmcs, nstep, mcs_start) # Получение данных из функции С++
All_Data.shape = ((int(((T_end-T_start)/T_step))+1,All_Data.size/(int(((T_end-T_start)/T_step))+1)))
np.savetxt('output/All_Data.dat', All_Data)

np.savetxt('output/E_mag.dat', np.transpose(np.array([All_Data[:,0], All_Data[:,1], All_Data[:,2]])))  #Кельвины (T, E1, E2)
np.savetxt('output/E_eV__J_mol).dat', np.transpose(np.array([All_Data[:,0], (All_Data[:,1]*kB)/Lattice.shape[0], (All_Data[:,1]*kB*96.32)/Lattice.shape[0]])))  # (T, E1 (meV), E2 (Дж/моль))
np.savetxt('output/M_x_y_z.dat', np.transpose(np.array([All_Data[:,0], All_Data[:,3], All_Data[:,5], All_Data[:,7], All_Data[:,4], All_Data[:,6], All_Data[:,8]])), fmt='%1.3f')  # (T, Mx1, My1, Mz1, Mx2, My2, Mz2)
np.savetxt('M_all.dat', np.transpose(np.array([All_Data[:,0], ((All_Data[:,3]**2 + All_Data[:,5]**2 + All_Data[:,7]**2)**0.5)])), fmt='%1.3f')  # (T, M_all)
'''

'''
spin_neighbors = np.zeros((J.shape[0], 3)) 
Neighbors = np.zeros((Lattice.shape[0], J.shape[0]))

for i in range(Lattice.shape[0]):
    for j in range(J.shape[0]):
        Neighbors[i, j] = -1 

for mri in range(Lattice.shape[0]):
    #print mri

    num_neighbors = 0
    curentspin = np.array([Lattice[mri,4], Lattice[mri,5], Lattice[mri,6]]) 

    IQ_curent_atom = Lattice[mri,0]
    x0 = Lattice[mri,1]
    y0 = Lattice[mri,2]
    z0 = Lattice[mri,3]

    for i in range(J.shape[0]):

        flag4 = False

        x = x0 + float(J[i,2])
        y = y0 + float(J[i,3])
        z = z0 + float(J[i,4])

        for j in range(Lattice.shape[0]):
            if (int(J[i,0]) == int(IQ_curent_atom) and int(J[i,1]) == int(Lattice[j,0])):# or (int(J[i,1]) == int(IQ_curent_atom) and int(J[i,0]) == int(Lattice[j,0])):
                if Lattice[j,1] == x and Lattice[j,2] == y and  Lattice[j,3] == z:
                    Neighbors[mri, num_neighbors] = j
                    flag4 = True
                    num_neighbors = num_neighbors + 1
        
       
        if flag4 == False:
            for k in range(26):
                for j in range(Lattice.shape[0]):
                    #x = x0 + L*primitive_vectors[k,0] + float(J[i,2])
                    #y = y0 + L*primitive_vectors[k,1] + float(J[i,3])
                    #z = z0 + L*primitive_vectors[k,2] + float(J[i,4])
                    
                    if k == 0:
                        x = x0 + L + float(J[i,2])
                        y = y0 + float(J[i,3])
                        z = z0 + float(J[i,4]) 
                    if k == 1:
                        x = x0 + float(J[i,2])
                        y = y0 + L + float(J[i,3])
                        z = z0 + float(J[i,4])
                    if k == 2:
                        x = x0 + float(J[i,2])
                        y = y0 + float(J[i,3])
                        z = z0 + L + float(J[i,4]) 
                    if k == 3:
                        x = x0 + L + float(J[i,2])
                        y = y0 + L + float(J[i,3])
                        z = z0 + float(J[i,4])       
                    if k == 4:
                        x = x0 + L + float(J[i,2])
                        y = y0 + float(J[i,3])
                        z = z0 + L + float(J[i,4])  
                    if k == 5:
                        x = x0 + float(J[i,2])
                        y = y0 + L + float(J[i,3])
                        z = z0 + L + float(J[i,4])  
                    if k == 6:
                        x = x0 + L + float(J[i,2])
                        y = y0 + L + float(J[i,3])
                        z = z0 + L + float(J[i,4])
                    if k == 7:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 + float(J[i,3]);
                        z = z0 + float(J[i,4]);
                    if k == 8:
    	  	        x = x0 + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 + float(J[i,4]);
                    if k == 9:
    	  	        x = x0 + float(J[i,2]);
                        y = y0 + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 10:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 + float(J[i,4]);
                    if k == 11:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 12:
    	  	        x = x0 + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 13:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);

                    if k == 14:
    	  	        x = x0 + L + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 + float(J[i,4]);
                    if k == 15:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 + L + float(J[i,3]);
                        z = z0 + float(J[i,4]);
                    if k == 16:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 + float(J[i,3]);
                        z = z0 + L + float(J[i,4]);
                    if k == 17:
    	  	        x = x0 + L + float(J[i,2]);
                        y = y0 + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 18:
    	  	        x = x0 + float(J[i,2]);
                        y = y0 + L + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 19:
    	  	        x = x0 + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 + L + float(J[i,4]);

                    if k == 20:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 + L + float(J[i,4]);
                    if k == 21:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 + L + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 22:
    	  	        x = x0 + L + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 23:
    	  	        x = x0 + L + float(J[i,2]);
                        y = y0 + L + float(J[i,3]);
                        z = z0 - L + float(J[i,4]);
                    if k == 24:
    	  	        x = x0 + L + float(J[i,2]);
                        y = y0 - L + float(J[i,3]);
                        z = z0 + L + float(J[i,4]);
                    if k == 25:
    	  	        x = x0 - L + float(J[i,2]);
                        y = y0 + L + float(J[i,3]);
                        z = z0 + L + float(J[i,4]);
                    

                    if (int(J[i,0]) == int(IQ_curent_atom) and int(J[i,1]) == int(Lattice[j,0])):# or (int(J[i,1]) == int(IQ_curent_atom) and int(J[i,0]) == int(Lattice[j,0])):
                        if Lattice[j,1] == x and Lattice[j,2] == y and  Lattice[j,3] == z:
                            #print mri, i, j
                            Neighbors[mri, num_neighbors] = j
                            num_neighbors = num_neighbors + 1
                    


np.savetxt('Neighbors.dat', Neighbors, fmt='%.0f')

colour = ['green', 'brown', 'violet', 'orange', 'violet']
size = [20, 20, 15, 15, 20]
fig=plt.figure()
ax = Axes3D(fig)

for i in range(Lattice.shape[0]):
    ax.plot([Lattice[i,1]], [Lattice[i,2]], [Lattice[i,3]], 'o', color = colour[int(Lattice[i,0])-1], markersize = size[int(Lattice[i,0])-1])

mri = 27
#print mri, Lattice[mri,:]

x0 = Lattice[mri,1]
y0 = Lattice[mri,2]
z0 = Lattice[mri,3]
ax.plot([x0,x0], [y0,y0], [z0,z0], 'or')

for j in Neighbors[mri,:]:
    if j != -1:
        x = Lattice[int(j), 1]
        y = Lattice[int(j), 2]
        z = Lattice[int(j), 3]

        ax.plot([x0,x], [y0,y], [z0,z], '-b', linewidth = 5)

ax.set_xlabel(u'X, м', size = 16)
ax.set_ylabel(u'Y, м', size = 16)
ax.set_zlabel(u'Z, м', size = 16)
ax.set_aspect(1)
plt.show()
'''
'''
mri = 0

colour = ['green', 'brown', 'violet', 'violet', 'violet']
size = [20, 10, 20, 20, 20]
fig=plt.figure()
ax = Axes3D(fig)

for i in range(Lattice.shape[0]):
    ax.plot([Lattice[i,1]], [Lattice[i,2]], [Lattice[i,3]], 'o', color = colour[int(Lattice[i,0])-1], markersize = size[int(Lattice[i,0])-1])

x0 = Lattice[mri,1]
y0 = Lattice[mri,2]
z0 = Lattice[mri,3]

ax.plot([x0,x0], [y0,y0], [z0,z0], 'or')

ax.plot([x0, x0+primitive_vectors[0,0]], [y0,y0+primitive_vectors[0,1]], [z0,z0+primitive_vectors[0,2]], '-r', linewidth = 5)
ax.plot([x0, x0+primitive_vectors[1,0]], [y0,y0+primitive_vectors[1,1]], [z0,z0+primitive_vectors[1,2]], '-g', linewidth = 5)
ax.plot([x0, x0+primitive_vectors[2,0]], [y0,y0+primitive_vectors[2,1]], [z0,z0+primitive_vectors[2,2]], '-b', linewidth = 5)

ax.set_xlabel(u'X, м', size = 16)
ax.set_ylabel(u'Y, м', size = 16)
ax.set_zlabel(u'Z, м', size = 16)
ax.set_aspect(1)
plt.show()
'''
