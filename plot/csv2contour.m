clear
clc
%close all

filename = '../N20_dt0p001_T0p01.csv';
Nx = 50;
Ny = 50;

%% READ AND SORT

A = csvread(filename,1,0);
A = sortrows(A,[1 2]);

%% RESHAPE TO MESHGRID FORMAT

X = reshape(A(:,1),Ny,Nx);
Y = reshape(A(:,2),Ny,Nx);
S = reshape(A(:,3),Ny,Nx);
V = reshape(A(:,4),Ny,Nx);

%% PLOT

figure();
contourf(X,Y,S);
title('Lid driven cavity problem - streamfuntion contour')
xlabel('x')
ylabel('y')
colorbar()

figure();
contourf(X,Y,V);
title('Lid driven cavity problem - vorticity contour')
xlabel('x')
ylabel('y')
colorbar()
