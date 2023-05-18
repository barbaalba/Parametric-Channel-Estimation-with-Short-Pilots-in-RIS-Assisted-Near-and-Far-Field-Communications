% This code is written to demonstrate how the beams stretched out of the
% rectangular border when elevation angle increases
close all; clear;clc;
initEl = [0,pi/8,pi/6]; % selected elevation angle
initAz = pi/2; % selected azimuth
M = 16; d = 1/4;
step = 1/(M*d);
figure('DefaultAxesFontSize',20,'defaultLineLineWidth',2,...
    'defaultAxesTickLabelInterpreter','latex');
for i = 1:length(initEl)
    x = (-3:step:3)/cos(initEl(i)) + sin(initAz);
    y = zeros(1,length(x))+ i/4;
    plot(x,y,'LineStyle','-','LineWidth',2,...
        'Marker','o','MarkerSize',10);
    hold on;
end
legend('$\theta=0$','$\theta=\pi/8$','$\theta=\pi/6$','interpreter','latex');
ylim([0,1]);
yticklabels([]);
grid on;
xlim([-4,4]);

% Figure Representation
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
