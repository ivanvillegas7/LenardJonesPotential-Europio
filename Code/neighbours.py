# -*- coding: utf-8 -*-
"""
Created on Mon Oct 18 15:04:05 2022

@author: Nicolás Gillard, Rubén Martínez and Iván Villegas

This program defines the function Neighbour_List(), which determines which are
the neighbouring atoms for each one of them. 
"""

from typing import List

import numpy as np

def Neighbour_List(N_atoms: int, Rl: float, Rc:float, N_FCC: int, Sigma: float,\
                   latcon: float, x: List[float], y: List[float], z: List[float]):

    '''
    This function requeries:
        ·The total number of atoms in the supercell (N_atoms)
        ·Cut off distance for the neighbour selection (Rl)
        ·Cut off distance for the potential (Rc)
        ·Number of FCC in the supercell (N_FCC)
        ·Deepness of the potential (Sigma)
        ·Lattice constant (lat_con)
        ·X-position of every atom in superlattice (x)
        ·Y-position of every atom in superlattice (y)
        ·Z-position of every atom in superlattice (z)

    With this arguments, the function saves which atom is neighbouring with any
    other atom (Neighbours). It also saves the index in which the indeces in which the
    neighbouring atoms start for the next atom (index).
    '''

    Rl = Rl/Sigma #In Sigma units

    POINT: List[int] = [] #Saves the position in the list

    LIST: List[int] = [] #The neighbours are stored here

    i_neighbours: int = 0 #Neighbours counter

    distance: List[float] = []

    BOXL: float = N_FCC * latcon/Sigma

    for i in range (N_atoms): #Loop over all the atoms

        atom_dist: List[float] = [] #Distances of each atom with its neighbours

        POINT.append(i_neighbours)

        for j in range(i+1, N_atoms):

            x_ij = x[i] - x[j] # Difference in X axis

            y_ij = y[i] - y[j] # Difference in Y axis

            z_ij = z[i] - z[j] # Difference in Z axis

            #Change the units to L=1

            x_ij = x_ij - round(x_ij/BOXL) * BOXL

            y_ij = y_ij - round(y_ij/BOXL) * BOXL

            z_ij = z_ij - round(z_ij/BOXL) * BOXL

            dist: float = np.sqrt(x_ij**2 + y_ij**2 + z_ij**2) #Compute de distance

            if dist <= Rl: #If it is a neighbour

                i_neighbours += 1

                LIST.append(i_neighbours)

                LIST[i_neighbours - 1]=j

                atom_dist.append(dist)

        distance.append(atom_dist)

    LIST: np.array = np.array(LIST)

    POINT[N_atoms - 1] = int(i_neighbours)
    
    '''
    #We can save the neighbour list of each atom at every instance to make sure
    #the programm works. We have disabled this option to optimise the program.
    
    #Save data in a txt file

    with open('neighbour_list.txt', 'w') as f:
        f.write(f'This is the file with the neighbour_list:\n')
        f.write(f'____________________________________________\n')
        f.write(f'Number of atoms: {N_atoms}\n')
        f.write(f'Cutoff distance: {2.7}\n')
        f.write(f'____________________________________________\n')
        for ind in range(0, N_atoms-1):
            f.write(f'Neighbour List of atom: {ind}\n')
            f.write(f'\n')
            f.write(f'Position of atom {ind}: {x[ind]}, {y[ind]}, {z[ind]}\n')
            f.write(f'Number of neighbours: {-(POINT[ind]-POINT[ind+1])}\n')
            f.write(f'\n')
            f.write(f'Index         Distances(Units L=1)\n')
            for i in range(POINT[ind], POINT[ind+1]):
                f.write(f'{LIST[i]}           {distance[ind][i-POINT[ind+1]]}\n')
    '''

    return POINT, LIST
