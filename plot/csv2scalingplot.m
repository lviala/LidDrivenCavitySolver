clear
clc
close all

filenames = {'../results/LDCoutput_Lx_1_Ly_1_Nx_100_Ny_100_Px_1_Py_1_T_50_Re_100.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_100_Ny_100_Px_2_Py_1_T_50_Re_100.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_100_Ny_100_Px_2_Py_2_T_50_Re_100.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_100_Ny_100_Px_4_Py_2_T_50_Re_100.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_100_Ny_100_Px_4_Py_4_T_50_Re_100.csv',...
            '../results/LDCoutput_Lx_1_Ly_1_Nx_100_Ny_100_Px_6_Py_5_T_50_Re_100.csv'};
Nx = 100;
Ny = 100;

LineStyleOrder = {'-k','--k','-.k',':k','-r','--r'};

perf = [];    

for ii = 1:length(filenames)
    %% READ AND SORT
    
    perf = [perf ; csvread(filenames{ii},1,0,[1 0 1 1])];
    
    A = csvread(filenames{ii},3,0);
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
    
    %% PLOT VELOCITY SLICES
    figure(1)
    hold on
    dy = (Y(3,1)-Y(2,1));
    idy = find((Y >= 0.5 - dy) & (Y < 0.5));
    plot(X(idy),velU(idy),LineStyleOrder{ii});
    
    figure(2)
    hold on
    dx = (X(1,3)-X(1,2));
    idx = find((X >= 0.5 - dx) & (X < 0.5));
    plot(Y(idx),velV(idx),LineStyleOrder{ii});
end

figure(3)
semilogx(perf(:,1),perf(1,2)./perf(:,2),'-ok')
xticks(perf(:,1))
title('Computational speedup compared to serial','interpreter','latex');
xlabel('Number of processors','interpreter','latex')
ylabel('Speedup','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
print(gcf,'ScalingPlot.png','-dpng','-r300');
