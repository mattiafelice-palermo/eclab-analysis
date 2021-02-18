#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:18:43 2021

@author: mattia
"""
import pandas as pd
import numpy as np
from scipy import integrate
import os
import sys

#%% Class for ECLab charge-discharge analysis


class CellCycling:
    """A class to read Biologic ECLab cell charge/discharge mpt files and
    compute relevant properties.
    
    ...
    
    Attributes
    ----------
    ncycles : int
        total number of charge/discharge cycles
    cretention : 1D numpy array of floats
        an array of the capacity retention of each charge/discharge cycle
    capacities : 1D numpy array of floats
        an array of the cell capacity for each charge/discharge cycle
    cyclelist : 1D numpy array of ints
        an array of the cycle indexes [1,2,...,ncycles] for use convenience
    efficiencies : 1D numpy array of floats
        an array of efficiency for each charge/discharge cycle
        
    Methods
    -------
    filter_row(column: str, value: float) -> pandas dataframe
        Returns a dataframe filtered by "value" in "column"
        
    """

    def __init__(self, folder, filename):
        """       

        Parameters
        ----------
        folder : str
            folder where the mpt file is contained
        filename : str
            name of the file to be read, including the mpt extension

        Returns
        -------
        None.

        """
        self._data = None  # dataframe with readings
        self._cycles = []  # list of Cycle objects
        self.columns = None # list of dataframe columns
        self.ncycles = None  # total number of cycles
        self.cretention = None  # list of capacity retentions for each cycle
        self.capacities = None  # list of capacities for each cycle
        self.cyclelist = None  # progressive index from 1 t ncycles for user
        self.efficiencies = {}  # dict of efficiencies for each cycle

        self.__read_mpt_file(folder, filename)  # read input file
        #self.__populate_cycles()  # create and add Cycle objects to _cycles[]
        self.__compute_cretention()  # compute cretention for each cycle
        self.__get_efficiencies()  # compute efficiency for each cycle

    def __read_mpt_file(self, folder, filename):
        path = os.path.join(folder, filename)

        with open(path, "r", encoding="utf8", errors="ignore") as f:
            
            delims=[] # contains cycle number, first and last line number
            beginning = None
            
            for line_num, line in enumerate(f):
                if "Number of loops : " in line:
                    self.ncycles = int(line.split(" ")[-1])
                
                
                # Before the output of the experiment, EClab lists the starting
                # and ending line of each loop. These will be used to slice
                # the pandas dataframe into the different cycles.
                if "Loop " in line:
                    loop_num = int(line.split(" ")[1])
                    first_pos = int(line.split(" ")[-3])
                    second_pos = int(line.split(" ")[-1])
                    delims.append([loop_num, first_pos, second_pos])
                    
                if "mode\t" in line:
                    beginning = line_num
                    break
            
            self._data = pd.read_table(
                    path, dtype=np.float64, delimiter = '\t', skiprows=beginning, decimal=","
            )
            
            self.columns = list(self._data.columns)
            self.cyclelist = np.arange(0, self.ncycles)
            
            cycle_num = 0
            
            # initiate Cycle object providing dataframe view within delims
            while cycle_num < self.ncycles:
                first_row = delims[cycle_num][1]
                last_row = delims[cycle_num][2]
                cycle = Cycle(cycle_num, self._data[first_row:last_row])
                self._cycles.append(cycle)
                cycle_num +=1

    def __getitem__(self, cycle):
        return self._cycles[cycle]

    def __iter__(self):
        for obj in self._cycles:
            yield obj

    def __compute_cretention(self, reference=1):
        # Retention is defined as the ration between the discharge for cycle=n
        # and cycle=1
        # The function also appends the cell capacities to the CellCycling obj
        initial_capacity = self._cycles[reference].capacity_discharge
        cretention = []
        capacities = []
        for cycle in self._cycles:
            cretention.append(cycle.capacity_discharge / initial_capacity)
            capacities.append(cycle.capacity_discharge)

        self.cretention = np.array(cretention)
        self.capacities = np.array(capacities)

    def __get_efficiencies(self):
        # Get efficiencies from Cycle objects and add it to CellCycling
        for L in ['C', 'V', 'E']:
            efficiency = [cycle.efficiencies[L] for cycle in self._cycles]
            self.efficiencies[L] = np.array(efficiency)

    def filter_row(self, column, value):
        """        
        
        Search for rows in dataframe that contain value in column
        and return the filtered dataframe.
        
        Parameters
        ----------
        column : str
            column of the dataframe that will be used for filtering
        value : float or int
            

        Returns
        -------
        Pandas Dataframe
            Pandas Dataframe containing the filtered rows.

        """
        return self._data[self._data[column] == value]


class Cycle:
    """
    Cycle is an object containing the cell readings for each cycle. The data
    is contained in a pandas dataframe with the following columns:
    '0. mode'
    'ox/red' # 1 during charge, 0 during discharge
    'error',
    'control changes',
     time/s' # timestep of the cell charge/discharge
    'control/V/mA',
    'Ewe/V' # Voltage
    'I/mA' # current
    'dq/mA.h',
    '(Q-Qo)/mA.h',
    'Q charge/discharge/mA.h',
    'Ece/V',
    'P/W',
     half cycle' # number of charge and discharge cycles
    'Q discharge/mA.h',
    'Q charge/mA.h',
    'Capacity/mA.h',
    'Efficiency/%',
    'control/V',
    'control/mA',
    'cycle number' # well, the cycle number
    'Ewe-Ece/V'
    
    The Cycle object also contains the efficiency and capacity (from discharge)
    of the cell cycle. 
    """

    def __init__(self, cycle, data):
        self._cyclen = cycle  # cycle number
        self._data = data  # view of CellCycling self._data pandas dataframe
        self.energy_charge = None
        self.energy_discharge = None
        self.capacity_charge = None
        self.capacity_discharge = None
        self.efficiencies = {}

        self.__compute_efficiencies()

    def __getitem__(self, key):
        """
        Operator overload. Allows access to dataframe columns using [] syntax.

        Parameters
        ----------
        key : str
            Pandas dataframe column name.

        Returns
        -------
        Pandas dataframe
            Returns the selected column of the cycle dataframe.

        """
        return self._data[key]

    def __compute_capacity(self):
        # Gets the last Capacity reading for the cycle discharge
        self.capacity_discharge = self.filter_row("ox/red", 0)["Capacity/mA.h"].iloc[-1]
        self.capacity_charge = self.filter_row("ox/red", 1)["Capacity/mA.h"].iloc[-1]
        
    def __compute_energy(self):
        energy = []
        for state in [0,1]: #0 is discharge, 1 is charge
            current = self.filter_row("ox/red", state)["(Q-Qo)/mA.h"]
            voltage = self.filter_row("ox/red", state)["Ewe/V"]        
            energy.append(integrate.trapezoid(voltage, x=current))
        
        self.energy_charge = energy[1]
        self.energy_discharge = -energy[0]
            

    def __compute_efficiencies(self):
        # Efficiency is computed as the ratio between the discharge and charge
        # energy. Energy is defined as C*V. Since C is constant for each
        # reading and simplifies in the ratio, efficienciy is computed as the
        # ratio between sum of voltages at discharge and charge.
        
        self.__compute_energy()
        self.__compute_capacity()
        
        coulomb_eff = self.capacity_discharge/self.capacity_charge
        energy_eff = self.energy_discharge/self.energy_charge
        self.efficiencies['C'] = coulomb_eff * 100
        self.efficiencies['E'] = energy_eff * 100
        self.efficiencies['V'] = energy_eff/coulomb_eff*100

    def filter_row(self, column, value):
        return self._data[self._data[column] == value]
