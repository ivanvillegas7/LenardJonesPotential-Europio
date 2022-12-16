# -*- coding: utf-8 -*-
"""

@author: Nicolás Gillard, Rubén Martínez and Iván Villegas

This program defines de function verlet_algorithm(), which uses the Verlet
integration to integrate the motion of the particles and calculate their
trajectories along each direction, following the Lennard Jones potential.
"""

import numpy as np

from neighbours import Neighbour_List as neig_l

from LennardJones import potentials

def verlet_algorithm(N_atoms, x, y, z, vx, vy, vz, Fx, Fy, Fz, D_Time, point,\
                     Neig_list, RC, RL, N_steps, Epsilon, Sigma, N_FCC, latcon):
    """
    This functions uses the Verlet integration to calculate the trajectories 
    of each particles. It computes the integration of the positions for each 
    time using the positions, velocities, and forces from the Lennard
    Jones potential.

    Parameters
    ----------
    N_atoms : int
        Total number of atoms.
    x : array
        Positions along the x axis at time t.
    y : array
        Positions along the y axis at time t.
    z : array
        Positions along the z axis at time t.
    vx : array
        Velocities along the x axis at time t.
    vy : array
        Velocities along the y axis at time t.
    vz : array
        Velocities along the z axis at time t.
    Fx : array
        Forces along the x axis at time t.
    Fy : array
        Forces along the y axis at time t.
    Fz : array
        Forces along the z axis at time t.
    D_Time: float
        Time step used in the integration Verlet algorithm.
    point: array
        List that points the position where the first neighbour of the atom
        can be found.
    Neig_list: array
        Neighbour list.
    RC : float
        Cut-off radius.
    Epsilon: float
        Variable of the Lennard-Jones potential.
    Sigma: float
        Variable of the Lennard-Jones potential.
    
    Returns
    -------
    positions_dt: array
        Position of the atom in the x, y and z axis at time t + D_Time.
    velocities_dt: array
        Velocities of the atom along the x, y and z axis at time t + D_Time
    point: array
        List that points the position where the first neighbour of the atom
        can be found.
    Neig_list: array
        Neighbour list.
    
    """
    
    #Initial acceleration
    ax = Fx
    ay = Fy
    az = Fz
    
    #Initialization for positions (t + D_Time)
    x_dt = np.zeros(N_atoms)
    y_dt = np.zeros(N_atoms)
    z_dt = np.zeros(N_atoms)
    
    #Initialization for velocities (t + D_Time/2)
    vx_dt_h = np.zeros(N_atoms)
    vy_dt_h = np.zeros(N_atoms)
    vz_dt_h = np.zeros(N_atoms)
    
    #Initialization for velocities (t + D_Time)
    vx_dt = np.zeros(N_atoms)
    vy_dt = np.zeros(N_atoms)
    vz_dt = np.zeros(N_atoms)
    
    #Initialization for dis(t + D_Time)
    
    for i in range(N_atoms):        
        
        #Compute the half step velocities
        vx_dt_h[i] = vx[i]+(1/2)*ax[i]*D_Time
        vy_dt_h[i] = vy[i]+(1/2)*ay[i]*D_Time
        vz_dt_h[i] = vz[i]+(1/2)*az[i]*D_Time
    
        #Propagate the positions from t to t + D_Time
        x_dt[i] = x[i] + vx_dt_h[i]*D_Time
        y_dt[i] = y[i] + vy_dt_h[i]*D_Time
        z_dt[i] = z[i] + vz_dt_h[i]*D_Time
        
        #Compute the forces from the positions at t + D_Time
        F_dt, E = potentials(N_atoms, point, Neig_list, x_dt, y_dt, z_dt, RC, \
                             Epsilon, Sigma, vx, vy, vz)
        
        fx_dt, fy_dt, fz_dt = F_dt[0], F_dt[1], F_dt[2]
        
        #Compute the acceleration (t+dt)
        ax_dt, ay_dt, az_dt = fx_dt/M, fy_dt/M, fz_dt/M
        
        #Compute the velocities(t+dt)
        vx_dt[i] = vx_dt_h[i] + ax_dt[i]*D_Time/2
        vy_dt[i] = vy_dt_h[i] + ay_dt[i]*D_Time/2
        vz_dt[i] = vz_dt_h[i] + az_dt[i]*D_Time/2
    
    '''
    #We can save the position of each atom at every instance to make sure the
    #programm works. We have disabled this option to optimise the program.
    
    #Print the positions at every instance
    
    with open('Positions_in_each_time_step.txt', 'w') as f1:
        f1.write(f'x                      y                         z\n')
        for i in range(108):
            f1.write(f'{x_dt[i]}    {y_dt[i]}    {z_dt[i]}\n')
    '''
    
    return [x_dt, y_dt, z_dt], [vx_dt, vy_dt, vz_dt], point, Neig_list
