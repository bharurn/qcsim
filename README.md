# QCSim
A quantum chemistry simulator. It is under initial stages of developemnt. Only single point energy calculation using Hartee-Fock method is available. 

The input should be given in the file 'data.integral'. It should contain the nuclear replusion energy, no. of atoms and no. of electrons in the system as the first line. Then the numerically pre-integrated values for the overlap integrals, core Hamiltonian and the two electron repulsion should be given in three blocks in that order, seperated by *. An example file for a water molecule with a bond length of 1.1 Å and a bond angle of 104°  with STO-3G basis set is given.

Integral calculation directly from basis set, frequency calculation, Hessian matrix calculation and geometry optimization will be added in due course to program.
