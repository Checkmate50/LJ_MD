## Lennard-Jones Molecular-Dynamics Simulator
Simple molecular dynamics simulator using Lennard-Jones potential.  Takes XYZ files as input and spits out XYZ files for use with VMD.

To compile, run `gcc lj_md.c`.  To run on Linux/Mac, use `./a.out simple.xyz simple`.  The two arguments to the code are the name of the input file and the format name of the output file (please don't include a file extension).  The results of running the code will be written to a series of XYZ files that represent a timestep animation of the atoms moving according to the LJ potential function.  These results can be viewed in the Visual Molecular Dynamics (VMD) simulator or a similar tool.
