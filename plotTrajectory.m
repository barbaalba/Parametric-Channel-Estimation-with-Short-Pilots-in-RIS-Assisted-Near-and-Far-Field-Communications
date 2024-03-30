function plotTrajectory(x_t,y_t,Az,El,conf,Xmax,Ymax,RIS_coor,d_bjo,d_NF)
%       RIS_Loc = [0,0,0]
%   ----------RIS---------- x = 0
%   |                      |
%   |                      |
%   |        ROOM          | 
%   |                      |
%   |                      |
%   ----------------------- x = Xmax
% -Ymax                   Ymax  
%
K = size(x_t,1); % number of user
T = size(x_t,2); % Time length
thet = 0:pi/180:2*pi;
xbjo = d_bjo*cos(thet) + RIS_coor(2);
ybjo = d_bjo*sin(thet) + RIS_coor(1);
xNF = d_NF*sin(thet) + RIS_coor(2);
yNF = d_NF*cos(thet) + RIS_coor(1);
for k = 1:K
    figure('defaultAxesFontSize',20,'DefaultLineLineWidth', 2,'defaultAxesTickLabelInterpreter','latex');
    fig = gcf;
    set(fig,'position',[60 50 1600 800]); % [left bottom width height]
    %tiledlayout(2,1);
    subplot(2,2,1);set(gca, 'ydir', 'reverse');
    xlabel('y','FontSize',20,'Interpreter','latex');ylabel('x','FontSize',20,'Interpreter','latex');
    grid on; xlim([-Ymax,Ymax]); ylim([0,Xmax]); hold on;
    plot(RIS_coor(2),RIS_coor(1),'Marker','square','MarkerSize',10,'Color','r'); % plot RIS in the room
    subplot(2,2,2);xlabel('Time (ms)','Interpreter','latex','FontSize',20);
    ylabel('Azimuth [degree]','Interpreter','latex','FontSize',20);grid on;ylim([-90,90]);
    xlim([1,T]);hold on;yticks([-90 -45 0 45 90]);
    subplot(2,2,[3,4]);xlabel('Time (s)','Interpreter','latex','FontSize',20);
    ylabel('Elevation [degree]','Interpreter','latex','FontSize',20); grid on;
    ylim([-90,0]);xlim([1,T]); hold on;
    for t = 2:T
    %nexttile(1);
    subplot(2,2,1);
    if strcmp(conf,'discrete')
        hold off;
        plot(RIS_coor(2),RIS_coor(1),'Marker','square','MarkerSize',10,'Color','r'); % RIS Location
        hold on;
        % The user movement
        plot(y_t(k,t-1:t),x_t(k,t-1:t),'Color','b','MarkerSize',20); 
        set(gca, 'ydir', 'reverse');
        xlabel('y','FontSize',20,'Interpreter','latex');
        ylabel('x','FontSize',20,'Interpreter','latex');
        xlim([-Ymax,Ymax]); ylim([0,Xmax]);grid on;
        plot(xbjo,ybjo,'--k','LineWidth',4);
        plot(xNF,yNF,'--k','LineWidth',4);
        hold off;
    else
        % the user trajectory
        plot(y_t(k,1:t),x_t(k,1:t),'Color','b');
        plot(xbjo,ybjo,'--k','LineWidth',4);
        plot(xNF,yNF,'--k','LineWidth',4);
    end
    %nexttile(2);
    subplot(2,2,2);
    plot(rad2deg(Az(k,1:t)),'r','Marker','*','MarkerSize',4);
    pause(0.00001);
%     subplot(2,2,[3,4]);
%     plot(rad2deg(El(k,1:t)),'k','Marker','*','MarkerSize',4);
%     pause(0.00001);
    end
end