# JJA-Simulation-Toolbox

Repository for MATLAB toolbox for simulation and analysis of Josephson junction arrays. Makes use of the JJAsim package by Martijn Lankhorst.

## Structure
* **functions** contains some functions that are used by the simulation and analysis scripts, such as generating array geometries and plotting critical current curves.
* **JJAsim** contains the JJAsim package that is responsible for the array calculations.
* **example scripts** contains base scripts for both simulation and analysis of critical current as a function of magnetic field for either a single array or a batch of different arrays, as well as a script for the visualization of vortex dynamics of an array. These serve as templates, and are to be copied and adjusted by the user.

## Basic usage
1. Make sure **functions** and **JJAsim** are added to the MATLAB path.
2. Copy one of the **example scripts** to your own folder, adapt it to your needs and run it.
