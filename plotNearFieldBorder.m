function plotNearFieldBorder(x_t,y_t,pltconf,Xmax,Ymax,RIS_coor,d_fraun)
K = size(x_t,1);
% circle border
thet = 0:pi/180:2*pi;
xunit = d_fraun*cos(thet) + RIS_coor(2);
yunit = d_fraun*sin(thet) + RIS_coor(1);
for k = 1:K
    figure('defaultAxesFontSize',20,'DefaultLineLineWidth', 2,'defaultAxesTickLabelInterpreter','latex');
    fig = gcf;
    set(fig,'position',[60 50 1600 800]); % [left bottom width height]
    set(gca, 'ydir', 'reverse');xlabel('y');ylabel('x');
    grid on; xlim([-Xmax,Xmax]); ylim([-Ymax,Ymax]); hold on; 
    plot(RIS_coor(2),RIS_coor(1),'Marker','square','MarkerSize',10,'Color','r');
    plot(y_t,x_t,'-b');
    plot(xunit,yunit,'--k','LineWidth',3);
end
end