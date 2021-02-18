# eclab-analysis
EClab-analysis is a Python library to read Biologic EClab mpt files.

##Installation
Download the EClab.py file and save it to your PYTHONPATH.

##Usage

```python
from EClab import CellCycling

folder='/yourpath'
file='filename.mpt'

readings=CellCycling(folder, file)

# If you want to know the number of cycles
print(readings.ncycles)
```

`readings` is a CellCyling object containing several `Cycle` objects. Individual
cycles can be accessed through square brackets operator `readings[x]`.
Each cycle object contains the raw data from the MPT file, returned as a numpy
array:

```python
first_cycle['time/s'] # returns the timesteps of each reading during the cycle
first_cycle['(Q-Qo)/mA.h'] # returns the current of each reading
first_cycle['Ewe/V'] # returns the voltage of each reading
```

The list of the keywords available for the raw data columns can be obtained
using `readings.columns`.

Each cycle object also contains some calculated properties:
```python
first_cycle.energy_charge
first_cycle.energy_discharge
first_cycle.capacity_charge
first_cycle.capacity_discharge
first_cycle.efficiencies['C'] # coulombic efficiency
first_cycle.efficiencies['V'] # voltage efficiency
first_cycle.efficiencies['E'] # energy efficiency
```

The `readings` object can be cycled in a foor loop e.g. to create the charge/discharge
capacity vs voltage plot:

```python
import matplotlib.pyplot as plt

for cycle in readings:
    plt.plot(cycle['(Q-Qo)/mA.h'], cycle['Ewe/V'])
    
# or if you wish to skip some of the cycles
for number in np.arange(0,9,2):
    plt.plot(readings[number]['(Q-Qo)/mA.h'], readings[number]['Ewe/V'])
```

Additionally, the `readings` object also conveniently contains numpy arrays
of the most useful calculated properties of the cycles, such as capacity retention,
energy, coulombic and voltage efficiencies. These can be used to plot their trend
as a function of the cycle number or time:

```python
plt.plot(readings.efficiencies['C'], label='Coulombic') # coulombic efficiency
plt.plot(readings.efficiencies['V'], label='Voltage') # voltage efficiency
plt.plot(readings.efficiencies['E'], label='Energy') # energy efficiency
plt.plot(readings.cretention[1:])
```
Note that the capacity retention is computed by default in reference to the
capacity of the second cycle.