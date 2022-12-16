# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 19:32:17 2022

@author: Nicolás Gillard Rubén Martínez and Iván Villegas

The following code is a function that reads an input.txt file regardless of
the provided units and returns them, so they are used in following codes.
"""

def txt_converter():
    """
    The file to be read MUST be called "input.txt". The way the .txt file must
    be given is the following: each line must be the variable followed by its value,
    leaving a space in between (for example: Length 10). It does not matter if they are 
    written with capital or lowercase letters or even if some abreviation sare used. 
    Default values will be asigned to the unprovided ones.
    
    -NumberFCCUnits or Number_FCC: int
        Number of FCC cells in each direction of the supercell of the simulation box
    -Reduced_density (adimensional): float
         Mass and sigma^ 3 divided by volume.
    -ReducedTemperature or Reduced_T (adimensional): float
        Reduced temperature inside the cube.
    -Sigma or s(in nm or m): float
        Deepness of the Lennard-Jones potential.        
    -Epsilon or e (in K or Cº): float = e/K_b
        Parameter of the Lennard-Jones potential.
    -CutOffPotential or RC (in sigma units): float
        CutOff distance for truncation of intermolecular potential. 
    -CutOffList or RL (in sigma units): float
        CutOff distance for the computation of Verlet neighbour list.
    -NumberOfSteps or Nsteps: int
        Number of steps in the molecular dynamic simulation, which, multiplied 
        by delta time gives us the total simulation time.
    -ReducedTimeStep or dt (in seconds): float
        Reduced time step of the molecular dynamic simulation.
    -WidthBinPairDistributionFunction or Width: float
        Width of the bins in the histogram to compute the pair distribution function.
    -Mass: float
        Atomic mass (in atomic mass units).

    """

#   Setup the default values for the different magnitudes and units

    Reduced_temperature_default         = 0.722
    Reduced_temperature                 = Reduced_temperature_default    

    ReducedTimeStep_default             = 0.001
    ReducedTimeStep                     = ReducedTimeStep_default

    Number_steps_default                = 100
    Number_steps                        = Number_steps_default  
    
    Number_FCC_default                  = 100
    Number_FCC                          = Number_FCC_default
    
    Epsilon_default                     = 1
    Epsilon_units_default               = "K"
    Epsilon                             = Epsilon_default 
    Epsilon_units                       = Epsilon_units_default  
    
    Sigma_default                       = 1
    Sigma                               = Sigma_default
    #Sigma_units                        = "nm"                               
        
    Reduced_density_default             = 0.8442
    Reduced_density                     = Reduced_density_default
    
    CutOffPotential_default             = 2.5*Sigma
    CutOffPotential                     = CutOffPotential_default  
    
    CutOffList_default                  = 2.7*Sigma
    CutOffList                          = CutOffList_default 
    
    Width_default                       = 0.1
    Width                               = Width_default
    
    Mass_default                        = 40
    Mass                                = Mass_default
    
    A_N_default                         = 18
    A_N                                 = A_N_default
    
    """
    What it will do is the following: 
    The "for" loop will read each line, it will then go through a series of "ifs" 
    and if the first word corresponds to a certain one in lowercase letters it will
    then go to the second word and get its value. If the value can be given in different
    units it will read the third word and change units if necessary.
    It has an "or" option if an abreviation is used.
    It is very similar to the one done in class with the ideal gas.
    """
    try:
        with open('input.txt', "r") as inputtxt:
            for line in inputtxt :
                # Lines are separated by blank spaces
                words = line.split()
                # If the line is only a blank line, it is discarded
                if len(words) > 1 :
                    #DELTA_TIME
                    if words[0].lower() == "reducedtimestep" or words[0].lower() == "dt":
                        ReducedTimeStep  = float(words[1])
                    
                   #NUMBER_STEPS
                    elif words[0].lower() == "numberofsteps" or words[0].lower() == "nsteps" :
                        Number_steps       = int(words[1])
                    
                    #NUMBER_FCC
                    elif words[0].lower() == "numberfccunits" or words[0].lower() == "number_fcc":
                        Number_FCC         = int(words[1])
                        
                    #EPSILON
                    elif words[0].lower() == "epsilon" or words[0].lower() == "e":
                        Epsilon            = float(words[1])
                        if len(words) > 2:
                            # Check the units of epsilon in the input file
                            Epsilon_units = words[2]
                            # If it is not in the SI system of units, transform it.
                            if Epsilon_units.upper() == "C" :
                                Epsilon = Epsilon + 273.15
                                Epsilon_units = Epsilon_units_default
                                
                    #SIGMA
                    elif words[0].lower() == "sigma" or words[0].lower() == "s":
                        Sigma              = float(words[1])
                        if len(words) > 2:
                          # Check the units of sigma
                          Sigma_units = words[2]
                          # If it is not in nm and it is in m, transform it.
                          if Sigma_units.upper() == "NM" :
                              Sigma = Sigma*1e-9
                        
                    #REDUCED_DENSITY
                    elif words[0].lower() == "reduceddensity":
                         Reduced_density = float(words[1])
                     
                    #CUTOFF_POTENTIAL
                    elif words[0].lower() == "cutoffpotential" or words[0].lower() == "rc":
                        CutOffPotential  = float(words[1])
                        
                    #CUTOFF_LIST
                    elif words[0].lower() == "cutofflist" or words[0].lower() == "rl":
                        CutOffList       = float(words[1])
                    
                    #WIDTH
                    elif words[0].lower() == "widthbinpairdistributionfunction" or words[0].lower() == "width":
                        Width            = float(words[1])
                    
                    #REDUCED_TEMPERATURE
                    elif words[0].lower() == "reducedtemperature" or words[0].lower() == "reduced_t":
                        Reduced_temperature = float(words[1])
                        
                    #MASS
                    elif words[0].lower() == "mass":
                        Mass = float(words[1])
                        
                    #ATOMIC NUMBER
                    elif words[0].lower() == "atomic_number" or words[0].lower() == "at_number" or words[0].lower() == "n_at":
                        A_N = float(words[1])

    except IOError:
        print("The file 'input.txt' does not exist")
    #In case you have misspelled the file name.
    
    #The number of atoms will then be:
    Number_of_atoms = 4 * (Number_FCC)**3

    #Print and return of the values 
    with open('output_conditions.txt', 'w') as f:
        f.write(f'Output file:\n')
        f.write(f'\n')
        f.write(f'Number of FCC units: {Number_FCC}; Number of FCC unit cells in each direction in the supercell\n')
        f.write(f'Which gives us a number of atoms :{Number_of_atoms} \n')
        f.write(f'Reduced density: {Reduced_density}; adimensional\n')
        f.write(f'Reduced temperature: {Reduced_temperature}; adimensional\n')
        f.write(f'Sigma: {Sigma} m\n')
        f.write(f'Epsilon: {Epsilon} K\n')
        f.write(f'Cutoff Potential: {CutOffPotential}; in units of sigma\n')
        f.write(f'Cutoff List: {CutOffList}; in units of sigma\n')
        f.write(f'Number of steps: {Number_steps} \n')
        f.write(f'Reduced time step: {ReducedTimeStep} s\n')
        f.write(f'Mass: {Mass} \n')
        f.write(f'Atomic number: {A_N} \n')
    
    return Sigma, Epsilon, Number_FCC, Reduced_density, CutOffPotential, \
           CutOffList, Number_steps, Reduced_temperature, Width, Number_of_atoms,\
           ReducedTimeStep, Mass, A_N
