# DirectStressMeasurementsEllipticBulging
This repository contains the files to directly determine the stresses during elliptic bulge tests, once the thin shell shape, loading and boundary conditions are known. The codes are developed in Matlab.

The folder contains two main files that work independently:

1) FDM_small_displacement_shell_fun.m - To directly calculate the stresses when the data is scattered, e.g. nodal data from FEA for an unstructured mesh

2) FDM_small_displacement_shell_gridded_fun.m - To directly calculate the stresses when the data is available in uniform 2D grids, e.g. data from digital image correlation (DIC) software

The remaining files are called by the above main functions. They are:
1) Principal_curvature_fun.m - called by FDM_small_displacement_shell_fun.m to compute principal curvatures and their directions

2) surfature.m - an open source code for calculating principal curvatures, called by Principal_curvature_fun.m
   ([Daniel Claxton (2022). Surface Curvature (https://www.mathworks.com/matlabcentral/fileexchange/11168-surface-curvature), MATLAB Central File Exchange. Retrieved August 18, 2022.]) 

3) derivate2.m - a function to calculate 2nd order numerical derivative for 1D arrays

4) derivative_2d.m - a function to calculate up to 4th order numerical derivatives for 2D arrays

Additionally, there are two example files that illustrate the application of the above functions:

1) Example_for_FEA.m - shows how FDM_small_displacement_shell_fun.m works

2) FDM_for_DIC.m - shows how FDM_small_displacement_shell_gridded_fun.m work

Users can open the above files independently and run them to see applications of the developed codes. The excel files contain data in gridded and scattered form, which are read by the abovementioned example files.


	
