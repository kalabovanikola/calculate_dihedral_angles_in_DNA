# calculate dihedral angles in DNA
## about
The program is used to calculate the angles of two valid nukleotides consecutively from given cif file. It can handle alt positions and heteroatoms of known purines and pyrimidines.
It can also reversely get pointers on the two sequence variants from given step ID. 

## usage
The main program takes the cif file name as an command line argument and returns a csv file with computed angles. 
The program for getting the step ID takes the step ID as the first argument and the name of the cif file as second and returns pointers on the sequence variants.
