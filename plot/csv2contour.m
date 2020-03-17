clear
clc
%close all

filename = '../results/LDCoutput_Lx_1_Ly_1_Nx_75_Ny_75_T_50_Re_100_conv1e7.csv';
Nx = 75;
Ny = 75;

%% READ AND SORT

A = csvread(filename,1,0);
A = sortrows(A,[1 2]);

%% RESHAPE TO MESHGRID FORMAT

X       = reshape(A(:,1),Ny,Nx);
Y       = reshape(A(:,2),Ny,Nx);
S       = reshape(A(:,3),Ny,Nx);
V       = reshape(A(:,4),Ny,Nx);
velU    = reshape(A(:,5),Ny,Nx);
velV    = reshape(A(:,6),Ny,Nx);

%% PLOT

figure();
subplot(2,2,1)
contourf(X,Y,S);
title('Lid driven cavity problem - streamfuntion contour')
xlabel('x')
ylabel('y')
colorbar()

subplot(2,2,2)
contourf(X,Y,V, linspace(-5,5,11));
title('Lid driven cavity problem - vorticity contour')
xlabel('x')
ylabel('y')
colorbar()

subplot(2,2,3)
contourf(X,Y,velU);
title('Lid driven cavity problem - U velocity')
xlabel('x')
ylabel('y')
colorbar()

subplot(2,2,4)
contourf(X,Y,velV);
title('Lid driven cavity problem - V velocity')
xlabel('x')
ylabel('y')
colorbar()
