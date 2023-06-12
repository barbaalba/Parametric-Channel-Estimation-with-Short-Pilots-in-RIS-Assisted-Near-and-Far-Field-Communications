close all; clear;clc;
%% this code is written to give us an intuition how red crosses are moved 
% in sine domain and fall outside the interval [-1,1] rectangle
M = 8; d = 1/2;
step = 1/(M*d); % sampling point in sine domain 
figure('DefaultAxesFontSize',20,'defaultLineLineWidth',2,'defaultAxesTickLabelInterpreter','latex');
hold on;
initEl = pi/6;
initAz = pi/2;
x = (-3:step:3)/cos(initEl) + sin(initAz);
y = zeros(1,length(x))+0.2*3;
plot(x,y,'LineStyle','-','LineWidth',2,'Marker','x',...
    'MarkerSize',10);
initAz = 0;
initEl = [pi/6 pi/8];
for i = 1:length(initEl)
    x = (-3:step:3)/cos(initEl(i)) + sin(initAz);
    y = zeros(1,length(x))+0.2*(length(initEl)-i+1);
    plot(x,y,'LineStyle','-','LineWidth',2,'Marker','x',...
        'MarkerSize',10);
end
x = (-3:step:3);
y = zeros(1,length(x));
plot(x,y,'LineStyle','-','LineWidth',2,'Color','b',...
    'Marker','o','MarkerSize',10);
grid on; box on;
xlim([-2,2]); 
ylim([-0.4,1]);
yticklabels([]);
legend('$\phi = \pi/6, \varphi = \pi/2$','$\phi=\pi/6, \varphi = 0$','$\phi=\pi/8, \varphi = 0$','$\phi = 0, \varphi = 0$','interpreter','latex');
% Figure Representation
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
