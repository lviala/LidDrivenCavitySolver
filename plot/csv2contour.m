clear
clc
close all

filename = '../results/LDCoutput_Lx_1_Ly_2_Nx_161_Ny_161_T_100_Re_100.csv';
Nx = 161;
Ny = 161;

%% READ AND SORT

A = csvread(filename,1,0);
A = sortrows(A,[1 2]);

%% STAT

expression = 'Re_(\d+).';
reynolds = regexp(filename,expression,'tokens');
[minS , id_minS] = min(A(:,3));
disp(['Re ', char(reynolds{1}), ...
    ': X = ' num2str(A(id_minS,1)),...
    ', Y = ' num2str(A(id_minS,2)),...
    ', S = ' num2str(A(id_minS,3)) ])

%% RESHAPE TO MESHGRID FORMAT

X       = reshape(A(:,1),Ny,Nx);
Y       = reshape(A(:,2),Ny,Nx);
S       = reshape(A(:,3),Ny,Nx);
V       = reshape(A(:,4),Ny,Nx);
velU    = reshape(A(:,5),Ny,Nx);
velV    = reshape(A(:,6),Ny,Nx);

%% PLOT

figure(1);
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

figure(2)
[C1,h1] = contour(X,Y,S,[linspace(-.10,0,11),-1e-2, -1e-4, 1e-7, 1e-6, 1e-5],'LineColor','k');
title('Streamfuntion Contour - Re = 100','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
clabel(C1,h1,-.10:.02:-.02,'interpreter','latex')
%print(gcf,'Stream_Ly2_Re100.png','-dpng','-r300');

figure(3)
[C2,h2] = contour(X,Y,V,linspace(-5,5,11),'LineColor','k');
t = clabel(C2,h2,'manual','interpreter','latex');
title('Vorticity Contour - Re = 100','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
xlabel('x','interpreter','latex')
ylabel('y','interpreter','latex')
%print(gcf,'Vort_Ly2_Re100.png','-dpng','-r300');

