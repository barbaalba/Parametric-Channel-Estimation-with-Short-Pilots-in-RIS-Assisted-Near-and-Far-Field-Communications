clear;clc;close all;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M = 32; d = 1/2;
W = Hcwgenerate(M,d,lambda);
fl = fieldnames(W);
resphi = -pi/2:pi/360:pi/2;
for k = 1:length(fl)
    w = W.(fl{k});
    gain =  abs(1/sqrt(M)*w'*ULA_Evaluate(lambda,M,resphi,d)).^2;
    figure;
    ax = polaraxes;
    polarplot(resphi+pi/2,pow2db(gain),'LineWidth',2); 
    hold on;
    polarplot(resphi+pi/2,repelem(0,1,length(resphi)),'LineStyle','--','Color','k','LineWidth',2);
    set(ax, 'ThetaTick', [0, 30, 60, 90, 120, 150, 180]);
    set(ax, 'TickLabelInterpreter', 'latex');
    set(ax, 'ThetaTickLabel',{'$-\frac{\pi}{2}$','$-\frac{\pi}{3}$','$-\frac{\pi}{6}$','$0$','$\frac{\pi}{6}$','$\frac{\pi}{3}$','$\frac{\pi}{2}$'});
    set(ax,'RLim',[-5 pow2db(M)]);
    set(gca,'fontsize',18);
    ax.ThetaLim = [0 180];
end




