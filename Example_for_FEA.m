clear;

% This code uses surfacture.m 
% [Daniel Claxton (2022). Surface Curvature (https://www.mathworks.com/matlabcentral/fileexchange/11168-surface-curvature), MATLAB Central File Exchange. Retrieved August 18, 2022.]

% Location of directory
path = pwd; 
filename = 'Coord_a_b_0.5_P_2.000_MPa_FEA.csv';

P = str2num(filename(:,17:21))*1e6;
 
a = 0.06;
b = 0.03;
v = 0.35; 
ho = 100e-6;

data  = xlsread( [path '\' filename],'L:O');  

x_inp=data(:,1);
y_inp=data(:,2);
z_inp=data(:,3);
h_inp=data(:,4);

[s1,s2,s12,w,R1,R2,x,y,z] = FDM_small_displacement_shell_fun(x_inp,y_inp,z_inp,h_inp,P,v,a,b,1);        
