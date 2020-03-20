clear
clc
close all

filenames = {'../results/LDCoutput_Lx_1_Ly_1_Nx_161_Ny_161_T_100_Re_100.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_161_Ny_161_T_100_Re_400.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_161_Ny_161_T_100_Re_1000.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_161_Ny_161_T_200_Re_3200.csv'};
Nx = 161;
Ny = 161;

LineStyleOrder = {'-k','--k','-.k',':k'};


for ii = 1:length(filenames)
    %% READ AND SORT
    
    A = csvread(filenames{ii},1,0);
    A = sortrows(A,[1 2]);
    
    %% STAT
    expression = 'Re_(\d+).';
    
    reynolds = regexp(filenames{ii},expression,'tokens');
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
    
    %% PLOT CONTOUR
    
    figure(ii);
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
    
    %% PLOT VELOCITY SLICES
    figure(length(filenames)+1)
    hold on
    idx = find(X == 0.5);
    plot(Y(idx),velU(idx),LineStyleOrder{ii});
    
    figure(length(filenames)+2)
    hold on
    idx = find(Y == 0.5);
    plot(X(idx),velV(idx),LineStyleOrder{ii});
end

figure(length(filenames)+1)
title('Horizontal velocity along the line x=0.5','interpreter','latex');
xlabel('y','interpreter','latex')
ylabel('u','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
legend('Re=100','Re=400','Re=1000','Re=3200','interpreter','latex','location','northwest')
print(gcf,'Uvel_x0p5.png','-dpng','-r300');

figure(length(filenames)+2)
title('Vertical velocity along the line y=0.5','interpreter','latex');
xlabel('x','interpreter','latex')
ylabel('v','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
legend('Re=100','Re=400','Re=1000','Re=3200','interpreter','latex')
print(gcf,'Vvel_y0p5.png','-dpng','-r300');
