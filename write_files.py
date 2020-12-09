# -*- coding: utf-8 -*-
#!/usr/bin/env python

import random
import os
import numpy as np

#### Сохранение файла задачи для кода C++ ####
def write_input(T_start, T_end, T_step, nmcs, nstep, mcs_start, L, sum_conc, num_atoms, J_shape, Lattice, J, magmom):

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
    file.write(str(sum_conc) + ' ')
    file.write(str(num_atoms)+ ' ')
    file.write(str(J_shape[0]) + ' ') #Число учитываемых обменных интеграллов
    file.close()

    file = open('input/Lattice.dat', "wb")
    np.savetxt(file, Lattice, fmt='%8f') #Данные о ячейке и спинах атомов
    file.close()
    
    file = open('input/J.dat', "wb")
    np.savetxt(file, J, fmt='%8f') 
    file.close()
    
    file = open('input/magmom.dat', "wb")
    np.savetxt(file, magmom, fmt='%8f') 
    file.close()
