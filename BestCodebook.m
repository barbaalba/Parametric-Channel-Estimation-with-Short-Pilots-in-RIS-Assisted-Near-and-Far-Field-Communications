%% This code is written to check if the codebook generated is orthogonal 
%% and plot a figure expressing what is the best initialization point and 
%% illustrate how the estimated dimension changes with initial angle pairs
clc;clear;close all;
%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 32; M_V = 32; M = M_H*M_V;
d_H = 1/2; d_V = 1/2; %In wavelengths
Th = -pi/2:pi/180:pi/2;
Az = -pi/2:pi/180:pi/2;
Num = zeros(length(Th),length(Az));
for thidx = 1:length(Th)
    parfor azidx = 1:length(Az)
        [ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,Az(azidx),Th(thidx));
        beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);
        x = abs(1/M*beamresponses'*beamresponses);
        orthnum = sum(x>0.5,'all');
%         if orthnum ~= CBL
%             break;
%         end
        Num(thidx,azidx) = CBL;
    end
end

figure('defaultAxesFontSize',20,'DefaultLineLineWidth', 2,...
    'defaultAxesTickLabelInterpreter','latex');
[X,Y] = meshgrid(Th,Az);
s = [-45 20];
sl = surfl(X,Y,Num',s,'light');
sl(2).Style = 'ambient';
xlabel('Elevation','FontSize',20,'Interpreter','latex');
sl(1).EdgeColor = 'none';
ylabel('Azimuth','FontSize',20,'Interpreter','latex');

% Find max angle pairs
mxval = max(Num,[],'all');
[phi_idx,varphi_idx] = find(Num == mxval);
linidx = find(Num == mxval);
Y = Th(phi_idx');
X = Az(varphi_idx');
hold on;
plot3(Y,X,Num(linidx),'LineStyle','none','Marker','o','MarkerSize',10,'Color','r',LineWidth=2);
yticks([-pi/2,-pi/4,0,pi/4,pi/2]);
xticks([-pi/2,-pi/4,0,pi/4,pi/2]);
xticklabels({'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'});
yticklabels({'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'});
