#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 11:18:43 2021

@author: mattia
"""
import pandas as pd
import numpy as np
from scipy import integrate

#%% Class for ECLab charge-discharge analysis


class GroupedCycles:
    """Class to analyze Biologic EClab output.

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

    def __init__(self, cycles):
        """

        Parameters
        ----------
        cycles : Cycle object
            list of Cycle objects

        Returns
        -------
        None.

        """
        self._cycles = cycles  # list of Cycle objects
        self.columns = None  # list of dataframe columns
        self.ncycles = None  # total number of cycles
        self.cretention = None  # list of capacity retentions for each cycle
        self.capacities = None  # list of capacities for each cycle
        self.cyclelist = None  # progressive index from 1 t ncycles for user
        self.efficiencies = {}  # dict of efficiencies for each cycle

        self.ncycles = len(self._cycles) # number of cycles

        self._remove_halfcycles() # delete cycles w/o (dis)charge capacity
        self.__compute_cretention()  # compute cretention for each cycle
        self.__get_efficiencies()  # compute efficiency for each cycle

    def __getitem__(self, cycle):
        return self._cycles[cycle]

    def __iter__(self):
        for obj in self._cycles:
            yield obj

    def _remove_halfcycles(self):
        for cycle in self._cycles:
            # charge or discharge values will be set to none if one of the two
            # is missing. It is enough to check for one
           
            if cycle.capacity_discharge is None:
                self._cycles.remove(cycle)
                #print("Removing cycle "+str(cycle._cyclen)+" missing (dis)charge")

        self.ncycles = len(self._cycles)

    def __compute_cretention(self, reference=1):
        # Retention is defined as the ration between the discharge for cycle=n
        # and cycle=1
        # The function also appends the cell capacities to the CellCycling obj
        try:
            initial_capacity = self._cycles[reference].capacity_discharge
        except:
            raise IndexError("Full cycles in MPT file: "+ str(self.ncycles)+". Cannot compute capacity retention.")
            
        cretention = []
        capacities = []
        for cycle in self._cycles:
            cretention.append(cycle.capacity_discharge / initial_capacity)
            capacities.append(cycle.capacity_discharge)

        self.cretention = np.array(cretention)
        self.capacities = np.array(capacities)

    def __get_efficiencies(self):
        # Get efficiencies from Cycle objects and add it to CellCycling
        for label in ["C", "V", "E"]:
            efficiency = [cycle.efficiencies[label] for cycle in self._cycles]
            self.efficiencies[label] = np.array(efficiency)

    def remove_incomplete_cycles(self):        
        # In case of partial charging followed by complete discharge, more
        # energy will be output than what was input. Therefore efficiency will
        # be higher than 100%. Remove these cycles and their corresponding
        # computed properties.
        for cycle in self._cycles:
            if cycle.efficiencies["E"] > 100:
                self._cycles.remove(cycle)
                #print("Removing cycle "+str(cycle._cyclen)+" due to eff > 100")
        self.__get_efficiencies()
        self.__compute_cretention()
        self.ncycles = len(self._cycles)


class CellCycling(GroupedCycles):
    '''
    Read multiple MPT files from Biologic EClab BCD routine.
    '''
    
    def __init__(self, *argv):
        '''
        

        Parameters
        ----------
        *argv : str
            List of filenames.

        Returns
        -------
        None.

        '''
        self._filenames = argv
        self._columnscheck = []

        all_cycles = []
        for filename in argv:
            all_cycles.extend(self.__read_single_mpt(filename))

        super().__init__(all_cycles)
        self.columns = self._columns_check()

    def __read_single_mpt(self, file):
        with open(file, "r", encoding="utf8", errors="ignore") as f:

            delims = []  # contains cycle number, first and last line number
            beginning = None

            ncycles = None
            for line_num, line in enumerate(f):
                if "Number of loops : " in line:
                    ncycles = int(line.split(" ")[-1])

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

            data = pd.read_table(
                file, dtype=np.float64, delimiter="\t", skiprows=beginning, decimal=","
            )

            self._columnscheck.append(list(data.columns))

            cycle_num = 0

            cycles = []

            # initiate Cycle object providing dataframe view within delims
            while cycle_num < ncycles:
                first_row = delims[cycle_num][1]
                last_row = delims[cycle_num][2]
                cycle = Cycle(cycle_num, file, data[first_row:last_row+1].copy())
                cycles.append(cycle)
                cycle_num += 1

            return cycles

    def _columns_check(self):
        columns = set.intersection(*[set(x) for x in self._columnscheck])
        return list(columns)


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

    def __init__(self, cycle, file, data):
        self._cyclen = cycle  # cycle number
        self._data = data  # view of CellCycling self._data pandas dataframe
        self.energy_charge = None
        self.energy_discharge = None
        self.capacity_charge = None
        self.capacity_discharge = None
        self.efficiencies = {}
        self.filename = file

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
        try:
            self.capacity_discharge = self.filter_row("ox/red", 0)[
                "Capacity/mA.h"
            ].iloc[-1]
            self.capacity_charge = self.filter_row("ox/red", 1)["Capacity/mA.h"].iloc[
                -1
            ]
            return True
        except:
            return False

    def __compute_energy(self):
        energy = []
        for state in [0, 1]:  # 0 is discharge, 1 is charge
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

        # check if cycle has done both charge and discharge, otherwise exit
        # without capacity we cannot compute coulombic efficiency
        if self.__compute_capacity() is False:
            return

        self.__compute_energy()

        coulomb_eff = self.capacity_discharge / self.capacity_charge
        energy_eff = self.energy_discharge / self.energy_charge
        self.efficiencies["C"] = coulomb_eff * 100
        self.efficiencies["E"] = energy_eff * 100
        self.efficiencies["V"] = energy_eff / coulomb_eff * 100

    def filter_row(self, column, value):
        return self._data[self._data[column] == value]
