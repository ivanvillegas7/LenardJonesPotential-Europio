# -*- coding: utf-8 -*-
"""
This file contains all the functions related with the initialization of the
atomic positions and atomic velocities.

Created on Sat Nov  19 16:30:46 2022

@author: Nicolás Gillard, Rubén Martínez and Iván Villegas
"""

import numpy as np

from typing import List

from random import Random

import matplotlib.pyplot as plt

def initial_position_FCC(Number_FCC: int, Reduced_temperature: float,\
                         Reduced_density: float, Sigma: float) -> List[np.array]:
    """
    This function generates the initial position of an atom for a supercell
    formed with N_C FCC cells.

    Parameters
    ----------
    Number_FCC: int
        Number of FCC units
    Reduced_temperature: float
        Reduced temperature.
    Reduced_density: float
        Reduced density.
    Sigma: float
        Sigma.
    Width : float
        Width of the histogram of initial velocities
    Returns
    -------
    Histogram of initial velocities.
    A txt file with the initial positions.
    A txt file with the initial velocities.
    Length : float
        Side of the unit cell.
    Positions : array
        Positions in the suppercell.
    latcon: float
        Lattice constat
    """

    #Lattice constant
    
    latcon=Sigma*(4/Reduced_density)**(1/3)
    
    #Positions in the FCC lattice
    
    x: List[float] = []
    
    y: List[float] = []
    
    z: List[float] = []
    
    Length: float = 1/Number_FCC
    
    cell_x=[0, 0.5*Length, 0, 0.5*Length]
    
    cell_y=[0, 0.5*Length, 0.5*Length, 0]
    
    cell_z=[0, 0, 0.5*Length, 0.5*Length]
    
    #Construct the lattice from the unit cell
    
    for iz in range(1, Number_FCC+1):
        
        for iy in range(1, Number_FCC+1):
            
            for ix in range(1, Number_FCC+1):
                
                for i in range(0, 4):
                    
                    x.append((cell_x[i]+Length*(ix-1)-0.5)*Number_FCC*\
                             latcon/Sigma)
                    
                    y.append((cell_y[i]+Length*(iy-1)-0.5)*Number_FCC*\
                             latcon/Sigma)
                    
                    z.append((cell_z[i]+Length*(iz-1)-0.5)*Number_FCC*\
                             latcon/Sigma)
                
    Position: List[np.array] =[x, y, z]
    
    #Print a file with the initial positions
    
    with open('output_initial_positions.txt', 'w') as f:
            f.write(f'This is a file with the initial positions:\n')
            f.write(f'__________________________________________\n')
        
            f.write(f'Index    Atomic coordinates(x, y, z)\n')
            for j in range(len(x)):
                f.write(f'{j}    {round(x[j],6)} {round(y[j],6)} {round(z[j],6)}\n')
    
    return Length, Position, latcon

def initial_velocities(Number_FCC: int, Reduced_temperature: float,\
                       Width: float, N_atoms: int) -> List[np.array]:
    """
    This function generates the initial velocities of an atom.
    The three cartesian coordinates are selected from a
    random distribution following the Maxwell-Boltzmann distribution.
    This corresponds with a normal distribution of mean value equal to zero,
    and standard deviation equal to
    sigma = sqrt((kB*T)/mass)
    Units: If mass in kg,
              Boltzmann constant in m^2 kg s^-2 K^-1
              and Temperature in K,
    then the velocities will be given in m/s.

    Parameters
    ----------
    Number_FCC : int
       Number of FCC lattices in the supercell
    Reduced_temperature  : float
       Reduced temperature
    N: int
       Number of atoms in the suppercell

    Returns
    -------
    initial_vel : numpy.array
       A numpy array of dimension 3, with the initial velocity of an atom
    """

    #Initialization of velocities
    
    v_x: List[float] = []
    
    v_y: List[float] = []
    
    v_z: List[float] = []
    
    var=np.sqrt(Reduced_temperature)
    
    g=Random(N_atoms)
    
    #Generate random velocities
    
    for _ in range(N_atoms):
        
        v_x.append(var*g.gauss(0, var))
        
        v_y.append(var*g.gauss(0, var))
        
        v_z.append(var*g.gauss(0, var))
        
    #Impose zero net moment
    
    nm_x=0
    
    nm_y=0
    
    nm_z=0
    
    for i in range(N_atoms):
        
        nm_x = nm_x+v_x[i]
        
        nm_y = nm_y+v_y[i]
        
        nm_z = nm_z+v_z[i]
        
    nm_x = nm_x/N_atoms
    
    nm_y = nm_y/N_atoms
    
    nm_z = nm_z/N_atoms
    
    for j in range(N_atoms):
        
        v_x[j]=v_x[j]-nm_x
        
        v_y[j]=v_y[j]-nm_y
        
        v_z[j]=v_z[j]-nm_z
       
    initial_vel: List[np.array] = [np.array(v_x), np.array(v_y), np.array(v_z)]
    
    #Chek: histogram of the velocities
    
    plt.figure()
    
    plt.hist(v_x, int(var/Width), label='Velocity axis X')
    
    plt.xlim(-2.5, 2.5)
    
    plt.legend()
    
    plt.savefig('output_velocity_X.pdf')
    
    plt.figure()
    
    plt.hist(v_y, int(var/Width), label='Velocity axis Y')
    
    plt.xlim(-2.5, 2.5)
    
    plt.legend()
    
    plt.savefig('output_velocity_Y.pdf')
    
    plt.figure()
    
    plt.hist(v_z, int(var/Width), label='Velocity axis Z')
    
    plt.xlim(-2.5, 2.5)
    
    plt.legend()
    
    plt.savefig('output_velocity_Z.pdf')
    
    #Print a file with the initial velocities
    
    with open('output_initial_velocities.txt', 'w') as f:
        f.write(f'This is a file with the initial velocities:\n')
        f.write(f'___________________________________________\n')
        
        f.write(f'Index    Velocities(v_x, v_y, v_z)\n')
        for j in range(N_atoms):
            f.write(f'{j}    {round(v_x[j],6)} {round(v_y[j],6)} {round(v_z[j],6)}\n')
    
    return initial_vel
