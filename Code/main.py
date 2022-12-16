# -*- coding: utf-8 -*-
"""
This is the main program of the melting solid.

Created on Wed Nov 23 15:55:36 2022

@author: Nicolás Gillard, Rubén Martínez and Iván Villegas
"""

import numpy as np

from typing import List

import matplotlib.pyplot as plt

"""
We import the functions needed, which are the initial positions and velocities
of the atoms, the conditions inside the box and the neighbour list.

Every output file or graph will be named 'output_{name}'.
"""

from inputs import txt_converter as inp

Sig: float

Eps: float

N_FCC: int

R_dens: float

r_c: float

r_l: float

N_steps: int

R_temp: float

Width: float

N: int

R_tStep: float
    
M: float

Sig, Eps, N_FCC, R_dens, r_c, r_l, N_steps, R_temp, Width, N, R_tStep, M = inp()

from init import initial_position_FCC as i_pos

L: float

Position: List[np.array]

latcon: float

from init import initial_velocities as i_vel

from neighbours import Neighbour_List as n_list

POINT: List[int]

LIST: List[int]

from LennardJones import potentials as pot

F: List[List[float]]

E_tot: float

Et: List[float] = []

from Verlet import verlet_algorithm as verlet

Pos_it: List[List[float]]

Vel_it: List[List[float]]

def main():

    L, Position, latcon = i_pos(N_FCC, R_temp, R_dens, Sig)

    x: np.array = Position[0]

    y: np.array = Position[1]

    z: np.array = Position[2]
    
    Velocities: List[np.array] = i_vel(N_FCC, R_temp, Width, N)

    v_x: np.array = Velocities[0]

    v_y: np.array = Velocities[1]

    v_z: np.array = Velocities[2]
    
    xt: List[List[float]] = [x]
    
    yt: List[List[float]] = [y]
    
    zt: List[List[float]] = [z]
    
    v_xt: List[List[float]] = [v_x]
    
    v_yt: List[List[float]] = [v_y]
    
    v_zt: List[List[float]] = [v_z]
    
    POINT, LIST = n_list(N, r_l, r_c, N_FCC, Sig, latcon, x, y, z)
    
    t: List[float] = np.linspace(0, N_steps*R_tStep, N_steps)
    
    movie = open("output_movie.AXSF", "w")
    
    movie.write(f'ANIMSTEPS     {N_steps}\n')
    
    movie.close()
    
    for i_step in range(N_steps):
        
        F, E_tot = pot(N, POINT, LIST, xt[i_step], yt[i_step], zt[i_step], r_c,
                       Eps, Sig, v_xt[i_step], v_yt[i_step], v_zt[i_step])
        
        F_x: List[float] = F[0]
        
        F_y: List[float] = F[1]
        
        F_z: List[float] = F[2]
        
        Et.append(float(E_tot))
        
        Pos_i, Vel_i, POINT, LIST = verlet(N, xt[i_step], yt[i_step], zt[i_step], 
                                           v_xt[i_step], v_yt[i_step], v_zt[i_step], 
                                           F_x, F_y, F_z, R_tStep, POINT, LIST,
                                           r_c, r_l, N_steps, Eps, Sig, N_FCC,
                                           latcon, M)
        
        xt.append(Pos_i[0])
        
        yt.append(Pos_i[1])
        
        zt.append(Pos_i[2])
        
        v_xt.append(Vel_i[0])
        
        v_yt.append(Vel_i[1])
        
        v_zt.append(Vel_i[2])        
        
        
        plt.figure()
        
        
        ax = plt.axes(projection='3d')
        
        ax.scatter3D(Pos_i[0], Pos_i[1], Pos_i[2], c='r', s = 50)
        
        ax.set_xlim(2,-3)
        
        ax.set_ylim(2,-3)
        
        ax.set_zlim(2,-3)
        
        ax.set_xlabel('x', labelpad=20)
        
        ax.set_ylabel('y', labelpad=20)
        
        ax.set_zlabel('z', labelpad=20)
        
        ax.set_title(f'Time step {i_step}')
        
        ax.grid()
        
        
        plt.savefig(f'Time step {i_step}.png')
        
        
        movie = open("output_movie.AXSF", "a")
        
        movie.write(f'ATOMS {i_step+1}\n')
        
        for i_atom in range(N):
            
            movie.write(f'1     {xt[i_step][i_atom]}     {yt[i_step][i_atom]}  {zt[i_step][i_atom]}\n')
            
        movie.close()
    
    plt.figure()
    
    plt.plot(t, Et, marker='o', label=r'Total energy of the system')
    
    plt.xlabel(r'$t$ [s]')
    
    plt.ylabel(r'$E$ [$\epsilon$]')
    
    plt.title(f'Conservation of energy\n')
    
    plt.legend()
    
    plt.grid(True)
    
    plt.savefig('output_energy_conservation.pdf')
    
    print('')
    
    print('The simulation is finished.')
    
    print('')
    
main()
