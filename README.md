# Force-Displacement-Analysis

Updated 13 September 2023
Author: Macy Mora-Antoinette

MATLAB R2019a  

<img width="400" alt="mts" src="https://github.com/macy-mora-antoinette/Force-Displacement-Analysis/assets/112992304/c535ed35-31ca-411a-acdd-b457682b31a1">

This matlab code is intended to help users compute the mechanical
properties from force-displacement data collected typically by three-point
bending of mouse femurs. 

It is separate into two parts:
PART 1: A simple for loop to iterate through all the raw data files
outputted by the MTS (mechanical testing system) located in Biomechanics
Labs (I believe the machine is owned by the van der Meulen group). It will
save the results in a csv.
PART 2: This function is called in PART 1 and is the main function that
computes the mechanical properties. It will save the plots of the
force-displacement data, which is good to have for the record. 

Key things to keep in mind:
1. This code is not automated. It requires user input for every data file.
The user will select the datapoints for stiffness, yeild force, etc. 
2. You must change the force to voltage equation in line 64. Typically,
a new calibration curve is created just before a round of testing, and you
try to use the same curve for the entire experiment by testing all your 
bones within a few days.
3. The columns for displacement and voltage can change depending on the
testing regime you pick on the MTS, so double check they're correct.
4. The FileName is the name you type in the MTS software when testing a new
specimen. This code requires a specific naming format. 
Ex C21-R-F-RF (mouseline-eartag-sex-rightfemur) all hyphens, NO SPACES
You can change this is you change the splitting mechanisms in lines 57 and 196
