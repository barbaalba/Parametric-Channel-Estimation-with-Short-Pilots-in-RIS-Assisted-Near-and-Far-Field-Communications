function plotTrajectory(x_t,y_t,Az,El,conf,Xmax,Ymax,RIS_coor)
K = size(x_t,1);
T = size(x_t,2);
for k = 1:K
    figure('defaultAxesFontSize',20,'DefaultLineLineWidth', 2,'defaultAxesTickLabelInterpreter','latex');
    fig = gcf;
    set(fig,'position',[60 50 1600 800]); % [left bottom width height]
    %tiledlayout(2,1);
    subplot(2,2,1);set(gca, 'ydir', 'reverse');xlabel('y');ylabel('x');
    grid on; xlim([-Xmax,Xmax]); ylim([-Ymax,Ymax]); hold on;
    plot(RIS_coor(2),RIS_coor(1),'Marker','square','MarkerSize',10,'Color','r');
    subplot(2,2,2);xlabel('Time (s)','Interpreter','latex','FontSize',20);
    ylabel('Azimuth','Interpreter','latex','FontSize',20);grid on;ylim([-90,90]);
    xlim([1,T]);hold on;yticks([-90 -45 0 45 90]);
    subplot(2,2,[3,4]);xlabel('Time (s)','Interpreter','latex','FontSize',20);
    ylabel('Elevation','Interpreter','latex','FontSize',20); grid on;
    ylim([-90,0]);xlim([1,T]); hold on;
    for t = 2:T
    %nexttile(1);
    subplot(2,2,1);
    if strcmp(conf,'discrete')
        hold off;
        plot(RIS_coor(2),RIS_coor(1),'Marker','square','MarkerSize',10,'Color','r'); % RIS Location
        hold on;
        % The user movement
        plot(y_t(k,t-1:t),x_t(k,t-1:t),'Color','b'); set(gca, 'ydir', 'reverse');xlabel('y');ylabel('x');xlim([-Xmax,Xmax]); ylim([-Ymax,Ymax]);grid on;
        hold off;
    else
        % the user trajectory
        plot(y_t(k,1:t),x_t(k,1:t),'Color','b');
    end
    %nexttile(2);
    subplot(2,2,2);
    plot(rad2deg(Az(k,1:t)),'r','Marker','*','MarkerSize',4);
    %pause(0.01);
    subplot(2,2,[3,4]);
    plot(rad2deg(El(k,1:t)),'k','Marker','*','MarkerSize',4);
    %pause(0.01);
    end
end