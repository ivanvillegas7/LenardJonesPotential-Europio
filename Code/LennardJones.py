# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 17:00:11 2022

@author: Nicolás Gillard, Rubén Martínez and Iván Villegas
"""

import numpy as np

from typing import List

def potentials(N, pointer, neighbourlist, x, y, z, r_c, Epsilon, Sigma, v_x, v_y, v_z):
    """
    This function uses the initial positions and velocities of the atoms, the conditions inside the 
    box and the neighbour list. By using the Lennard-Jones potential, it will obtain the forces 
    of each atom and the kinetic, potential and total energy.
    
    Parameters- which are the same as shown in the other functions really-.
    ----------
    N : int
        Total number of atoms.
    pointer (also referred as POINT) : array
        List that points the position where the first neighbour of the atom can be found.
    neighbourlist (also referred as LIST) : array
        Neighbour list.
    -Rc: float
        Cutoff potential: limit distance, a larger distances potential value will be taken as 0. 
    -Sigma: float
        Deepness of the potential
    -Epsilon: float = e/K_b
        Energy value to make everything dimensionless
    -x, y, z : array
        Positions along the x, y and z axis, respectively.
    -v_x, v_y, v_z: array
        Velocities along the x, y and z axis, respectively.
    Returns
    -------
    Three txt files, with the potential, kinetic and total energy.
    V : array
        Lennard-Jones potential.
    K : array
        The kinetic energy.
    E : array
        The total energy, which is just the sum of both previous ones.
    F : array
        The force of each atom along the x, y and z axis, respectively.
    """
        
    #Total distance
    r: List[float] = []
    
    for i in range(N):
        r.append(np.sqrt(x[i]**2+y[i]**2+z[i]**2))
    
    """
    We know create the empty vectors which will be filled with the Lennard-Jones
    potential and the forces. We create one force vector for each dimension, and 
    we will join them after.
    """
    
    V   = np.zeros(N) #Lennard Jones potential
    F_x = np.zeros(N) #Force vector on the x axis
    F_y = np.zeros(N)
    F_z = np.zeros(N)
    
    #Loop over all atoms except the last one
    for i in range(N-1):
        
        #Loop from the first to the last neighbour of atom i
        for j in range(int(pointer[i]), int(pointer[i+1]-1)):
            
            #Relative position between atom i and j 
            posi = r[i]
            posj = r[neighbourlist[j]]
            relative_pos = round(posj-posi, 10)
            
            #If this relative position is negative, we take its absolute value.
            if relative_pos>0:
                distance =np.sqrt(relative_pos)
                
            else:
                distance=-np.sqrt(abs(relative_pos))
            
            '''
            We now want to neglect atoms which are further away from the cut-off
            radius, where the potential (and so the forces) will be taken as 0.
            '''
  
            if distance > r_c or distance == 0:
                V[i]   = 0
                F_x[i] = 0
                F_y[i] = 0
                F_z[i] = 0
                
            else:
                #Potential
                Vi   = 4*Epsilon*((Sigma/distance)**12-(Sigma/distance)**6)
                V[i] = V[i]+Vi
                
                #Force
                f = (24*Epsilon/distance)*(2*(Sigma/distance)**12-(Sigma/distance)**6)
                
                #X-Force 
                fxi    = f*(x[neighbourlist[j]]-x[i])
                F_x[i] = F_x[i]+fxi
                
                #Y-Force
                fyi    = f*(y[neighbourlist[j]]-y[i])
                F_y[i] = F_y[i]+fyi
                
                #Z-Force 
                fzi    = f*(z[neighbourlist[j]]-z[i])
                F_z[i] = F_z[i]+fzi
    
    #We now create a force vector which includes the three dimensions
    F : List[np.array] = [F_x, F_y, F_z]
    
    #Kinetic energy, which is calculated just be using the results provided by init.py
    #and total energy, which is the sum of both energies for each value.
    velocity = np.zeros(N)
    K        = np.zeros(N)
    E        = np.zeros(N)
    E_tot: float = 0
    
    for i in range(N):
        velocity[i] = np.sqrt(v_x[i]**2+v_y[i]**2+v_z[i]**2)
        K[i]        = (velocity[i]**2)/2
        #Total energy
        E[i]        = V[i] + K[i]
        E_tot           = E_tot + E[i]
    
    #We save the potential, kinetic and total energy in .txt files
    with open('potentials.txt', 'w') as pot:
        
        pot.write(f'Atom    Potential energy\n')
        
        for i in range(N):
            
            pot.write(f'{i+1}         {V[i]}\n')
        
    with open('kinetic_energy.txt', 'w') as kin:
        
        kin.write(f'Atom    Kinetic energy\n')
        
        for i in range(N):
            
            kin.write(f'{i+1}         {K[i]}\n')
            
    with open('total_energy.txt', 'w') as tot:
        
        tot.write(f'Atom    Total energy\n')
        
        for i in range(N):
            
            tot.write(f'{i+1}         {E[i]}\n')
    
    return F, E_tot
