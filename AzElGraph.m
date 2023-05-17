%% AZ-EL Graph 
% this code is written to produce a plot of all the angles that results in
% orthogonal beams
M_V = 8; M_H = 8; d = 1/4;
[theta,phi,num] = UPA_BasisElupnew(M_V,M_H,d,d,pi/2,0);
% plot the results
figure('DefaultAxesFontSize',20,'defaultLineLineWidth',2,'defaultAxesTickLabelInterpreter','latex')
x = linspace(-90,90,360);
fn = fieldnames(phi);
for i = 1:length(theta)
    y = theta(i);
    plot(x,repelem(rad2deg(y),length(x)),'b'); % blue lines
    hold on;
    for j = 1:length(phi.(fn{i}))
        plot(rad2deg(phi.(fn{i})(j)),rad2deg(y),...
            'r','Marker','x','MarkerSize',15); % red cross
    end
end
% configure the figure
xlim([-90,90]);
ylim([-90,90]);
xlabel('Azimuth','Interpreter','latex');
ylabel('Elevation','Interpreter','latex');
xticks(-90:20:90);
xticklabels({'$-90$','$-70$','$-50$','$-30$','$-10$','$10$','$30$','$50$','$70$','$90$'});
yticks(-90:20:90);
yticklabels({'$-90$','$-70$','$-50$','$-30$','$-10$','$10$','$30$','$50$','$70$','$90$'});
grid on;
% Figure size
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
