#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:18:43 2021

@author: mattia
"""
import pandas as pd
import numpy as np
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
        self.ncycles = None  # total number of cycles
        self.cretention = None  # list of capacity retentions for each cycle
        self.capacities = None  # list of capacities for each cycle
        self.cyclelist = None  # progressive index from 1 t ncycles for user
        self.efficiencies = None  # list of efficiencies for each cycle

        self.__read_mpt_file(folder, filename)  # read input file
        self.__populate_cycles()  # create and add Cycle objects to _cycles[]
        self.__compute_cretention()  # compute cretention for each cycle
        self.__get_efficiencies()  # compute efficiency for each cycle

    def __read_mpt_file(self, folder, filename):
        path = os.path.join(folder, filename)

        with open(path, "r", encoding="utf8", errors="ignore") as f:
            for line in f:
                if "Number of loops : " in line:
                    self.ncycles = int(line.split(" ")[-1])
                    break

            self._data = pd.read_table(
                path, dtype=np.float64, delimiter="\t", skiprows=410, decimal=","
            )

    def __populate_cycles(self):
        # Filter CellCycling  pandas dataframe self._data by cycle number and
        # creates an Cycle object containing readings only for that cycle
        cycle = 0
        cyclelist = []
        while cycle < self.ncycles:
            data = self._data[self._data["cycle number"] == cycle]
            self._cycles.append(Cycle(cycle, data))
            cyclelist.append(cycle)
            cycle += 1

        self.cyclelist = np.array(cyclelist)

    def __getitem__(self, cycle):
        return self._cycles[cycle]

    def __iter__(self):
        for obj in self._cycles:
            yield obj

    def __compute_cretention(self):
        # Retention is defined as the ration between the discharge for cycle=n
        # and cycle=1
        # The function also appends the cell capacities to the CellCycling obj
        initial_capacity = self._cycles[0].capacity
        cretention = []
        capacities = []
        for cycle in self._cycles:
            cretention.append(cycle.capacity / initial_capacity)
            capacities.append(cycle.capacity)

        self.cretention = np.array(cretention)
        self.capacities = np.array(capacities)

    def __get_efficiencies(self):
        # Get efficiencies from Cycle objects and add it to CellCycling
        efficiencies = [cycle.efficiency for cycle in self._cycles]
        self.efficiencies = np.array(efficiencies)

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
    Cycle is and object containing the cell readings for each cycle. The data
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
        self._data = data  # pandas dataframe
        self.efficiency = 0
        self.capacity = 0

        self.__compute_efficiency()
        self.__get_capacity()

    def __getitem__(self, key):
        return self._data[key]

    def __compute_efficiency(self):
        # Efficiency is computed as the ration between the discharge and charge
        # energy. Energy is defined as C*V. Since C is constant for each
        # reading and simplifies in the ratio, efficienciy is computed as the
        # ratio between sum of voltages at discharge and charge.
        Echarge = self.filter_row("ox/red", 1)["Ewe/V"].sum()
        Edischarge = self.filter_row("ox/red", 0)["Ewe/V"].sum()
        self.efficiency = Edischarge / Echarge * 100

    def __get_capacity(self):
        # Gets the last Capacity reading for the cycle discharge
        self.capacity = self.filter_row("ox/red", 0)["Capacity/mA.h"].iloc[-1]

    def filter_row(self, column, value):
        return self._data[self._data[column] == value]
