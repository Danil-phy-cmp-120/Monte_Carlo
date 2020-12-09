# -*- coding: utf-8 -*-
#!/usr/bin/env python

import random
import os
import numpy as np

import read_files
import write_files

'''
L = 10
d_a = 2.0

T_start = 10.0
T_end = 1600.0
T_step = 20.0

nmcs = 100000
nstep = 400
mcs_start = int(nmcs/10)
'''

L = input('Enter size of lattice for calculation (L x L x L):\n')
L = int(L)

d_a  = input('Enter the maximum distance to take into account the magnetic exchange integrals:\n')
d_a = float(d_a)

T_start = input('Enter the minimum temperature:\n')
T_end = input('Enter the maximum temperature:\n')
T_step = input('Enter the temperature step:\n')
T_start = float(T_start); T_end = float(T_end); T_step = float(T_step)

nmcs = input('Enter the number of steps Monte Carlo:\n')
nstep = input('Enter the number of steps through which averaging:\n')
nmcs = int(nmcs); nstep = int(nstep)
mcs_start = int(nmcs/10)


primitive_vectors, A, basis = read_files.read_cell()
IQ, conc = read_files.read_atoms()
J = read_files.read_J(d_a)

for i in range(J.shape[0]): # Перевод векторов из Декатового базиса в базис решетки
    J[i,2:5] = A.dot(J[i,2:5])

for i in range(len(conc*L**3)): # Небольшое изменение концентрации под текущую ячейку чтобы число атомов было int
    conc[i] = round(conc[i]*L**3)/(L**3)

#### Премешивание атомов на одной позиции ####
atom_index = 1
atoms = []
atoms_cpa = []
for j in range(len(IQ)):
    if conc[j] != 1.0:
        for m in range(int(L**3*conc[j])):
            atoms_cpa += [atom_index]
    if conc[j] == 1.0 or j == len(IQ)-1:
        random.shuffle(atoms_cpa) 
        atoms = atoms + atoms_cpa
        atoms_cpa = []
        for m in range(L**3):
            atoms += [atom_index]

    if j < len(IQ)-1:
        if IQ[j] != IQ[j+1]:
            atom_index = atom_index + 1

#### Транслируем атомы базиса по векторам трансляций ####     
Lattice = np.zeros((int(sum(conc))*L**3, 7))               
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
                Lattice[n,6] = 1 
                n = n + 1

Lattice = sorted(Lattice, key=lambda a_entry: a_entry[0])
Lattice = np.array(Lattice)

for i in range(Lattice.shape[0]): # Перевод векторов из Декатового базиса в базис решетки
    Lattice[i,1:4] = A.dot(Lattice[i,1:4])


magmom = read_files.read_magmom(max(atoms))

write_files.write_input(T_start, T_end, T_step, nmcs, nstep, mcs_start, L, int(sum(conc)), max(atoms), J.shape, Lattice, J, magmom)
