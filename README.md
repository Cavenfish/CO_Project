# CO Project

This repository contains all files used in the CO ice energy
dissapation research project.

## CO Project V 0.0

This contains a tarball of the orginal source code and scripts
written by Fabian Sies. 

## AutoMD

This is a python package that improves workflow while running simulations 
related to the CO Project. The config file contains several helpful
functions that can be used to automate the process of running MD calculations

### prep_system

**`prep_system(xyz)`**

```
This function turns and xyz file into an Atoms object with MvH calculator.

Parameters
----------
xyz: string
	A string pointing to an xyz file location.

Returns
-------
system: ASE Atoms Object
	A object representing the xyz system given to the function.
calc: ASE Calculator Object
	This will be always be the MvH calculator object with the system
	assigned to it.

Examples
--------
>>> import AutoMD as md
>>>
>>> xyz = './you_file.xyz'
>>> 
>>> system, calc = md.prep_system(xyz)
>>>
```

### stretch_molecule

**`stretch_molecule(xyz, swap, masses, r)`**

```
Stretches the bond length (without moving center of mass) of a molecule, to
give it a desired amount of vibrational energy. ONLY WORKS FOR DIATOMIC

Parameters
----------
xyz: string
        A string pointing to an xyz file location.
swap: list 
	Indices of the atoms being moved
masses: list
	Masses of the atoms being moved
r: float
	Desired bond length after stretch

Notes
-----
'swap' and 'masses' must have the same length and their indicies should 
corrspond to the same atoms.

This function will automatically save a new xyz file with the excited 
molecule within the system given.

Function will always save new file with '_excited' appended to file name

Returns
-------
new_name: string
	A string pointing to an xyz file of the system with the
	newly excited molecule within it.

Examples
--------
>>> import AutoMD as md
>>>
>>> xyz    = './your_file.xyz'
>>> swap   = [ 0,  1]
>>> masses = [12, 16]
>>> r      = 1.212
>>>
>>> new_name = md.stretch_molecule(xyz, swap, masses, r)
>>>
>>> new_name
'./your_file_excited.xyz'
```

### Morse_excitation

**`Morse_excitation(nu, n)`**

```
Solves Morse potential equation to determine the bond length of an excited
molecule. PARAMETERS CURRENTLY ONLY SUPPORT CO

Parameters
----------
nu: float
	Vibrational frequency of molecule being excited (in wavenumbers cm^-1)
n: integer
	Excitation number 

Returns
-------
r: float
	Bond length of excited molecule

Examples
--------
>>> import AutoMD as md
>>> 
>>> nu = 2100
>>> n  = 1
>>>
>>> r  = md.Morse_excitaton(2100, 1)
>>>
>>> r
1.212
```

### geo_opt

**`geo_opt(xyz)`**

```
Runs a BFGS geometry optimization on the system given

Parameters
----------
xyz: string
        A string pointing to an xyz file location.

Returns
-------
opt_f: string
	A string pointing to the optimized xyz file location

Examples
--------
>>> import AutoMD as md
>>>
>>> xyz   = './your_file.xyz'
>>>
>>> opt_f = md.geo_opt(xyz)
>>>
>>> opt_f
'./your_file_opt.xyz'
```

### calc_vibs

**`calc_vibs(xyz)`**

```
Calculates the vibrational frequencies of a given system

Parameters
----------
xyz: string
        A string pointing to an xyz file location.

Returns
-------
None

Examples
--------
>>> import AutoMD as md
>>> 
>>> xyz = './your_file.xyz'
>>>
>>> md.calc_vibs(xyz)
```

### add_isotope

**`add_isotope(xyz, pos, masses)`**

```
Adds a given isotope to a given system

Parameters
----------
xyz: string
        A string pointing to an xyz file location.
pos: list
	Indices of atoms being replaced
masses: list
	Masses of the isotopic atoms

Returns
-------
new_name: string
	A string pointing to the newly made xyz file

Examples
--------
>>> import AutoMD as md
>>>
>>> xyz      = './your_file.xyz'
>>> pos      = [ 0,  1]
>>> masses   = [12, 18]
>>>
>>> new_name = md.add_isotope(xyz, pos, masses)
>>>
>>> new_name
'./your_file_isotope.xyz'  
```



## Source Code

This contains the source code for the ASE calculator.

## XYZ Files
