# -*- coding: utf-8 -*-
#!/usr/bin/env python

import random
import os
import numpy as np

#### Считывание файла JXC.out ####
def read_cell(path='*JXC.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()

    #### Считываем вектора трансляций, базис ####
    Data_lattice = []
    flag = False
    for line in lines:
        inp = line.split()
        if len(inp) == 0:
            continue
        if flag == True and len(inp) == 5 and inp[0] == '(':
           Data_lattice += [float(inp[1].replace(',','')), float(inp[2].replace(',','')), float(inp[3])]
        if inp[0] == '<INIT_MOD_LATTICE>' and len(inp) >= 1:
            flag = True
        if inp[len(inp)-1] == '2*pi/a' and len(inp) > 1:
            break

    primitive_vectors = np.array(Data_lattice[:9])
    primitive_vectors.shape = (3, 3)

    A = np.linalg.inv(primitive_vectors) # Матрица перехода от декартового базиса к текущему
 
    basis = np.array(Data_lattice[9:len(Data_lattice)])
    basis.shape = (int(len(basis)/3), 3)

    return primitive_vectors, A, basis

#### Считываем тип атомов ####
def read_atoms(path='*JXC.out'):
    f = open(path, "r")
    lines = f.readlines()
    f.close()
    
    IQ = []
    conc = []
    flag = False
    for line in lines:
        inp = line.split()
        if flag == True  and len(inp) >= 9:
            for i in range(int(inp[len(inp)-1]) - int(inp[8]) + 1):
                IQ += [inp[1]]
                conc += [float(inp[7])]
        if len(inp) > 3:
            if inp[0] == 'type' and inp[1] == 'TXTT' and inp[2] == 'NL':
                flag = True
        if len(inp) == 1:
            if inp[0] == 'mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm':
                break
    conc = np.array(conc)
    return IQ, conc


#### Считываем обменные интегралы ####
def read_J(d_a, path='*JXC.out'):

    f = open(path, "r")
    lines = f.readlines()
    f.close()

    J = []
    flag = False
    for line in lines:
        inp = line.split()
        if len(inp) == 0:
            continue

        if flag == True and len(inp) == 11 and abs(float(inp[10])*1000) > 0.01 and float(inp[8]) <= d_a:
            J += [IQ, IT, float(inp[5]), float(inp[6]), float(inp[7]), float(inp[8]), float(inp[10])*1000]
            J += [IT, IQ, float(inp[5]), float(inp[6]), float(inp[7]), float(inp[8]), float(inp[10])*1000]
            flag = False

        if inp[0] == 'ITAUIJ' and inp[1] == 'ITAUJI' and len(inp) >= 1:
            flag = True

        if inp[0] == 'IQ' and inp[6] == 'JQ' and len(inp) == 12 and flag == False:
            IQ = float(inp[5]); IT = float(inp[11])

    J = np.array(J)
    J.shape = (int(J.size/7), 7)
    J = np.array(list(map(list, {tuple(x) for x in J})))
    J = sorted(J, key=lambda x: (x[5], -x[6]))
    J = np.array(J)

    return J


#### Считывание магнитных моментов из SCF.out ####
def read_magmom(num_atoms, path='*SCF.out'):

    f = open(path, "r")
    lines = f.readlines()
    f.close()

    magmom = np.zeros(num_atoms)
    flag = False
    n = 0
    for i in range(num_atoms):
        for line in lines[len(lines)-200:len(lines)]:
            inp = line.split()
            if len(inp) > 5 and inp[1] == 'E=' and inp[4] == 'IT=' and int(inp[5]) == i+1:
                flag = True
            if len(inp) > 9 and flag == True and inp[0] == 'sum':
                magmom[n] = float(inp[4])
                flag = False
        n = n + 1
    return magmom
