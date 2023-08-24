% This code is written to illustrate the phase variations over angtennas
clear;clc;close all;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_V = 512; M_H = 512; M=M_H*M_V; d_V = 1/32; d_H = 1/32;
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D = sqrt(Hsize^2+Vsize^2); % Diagonal size of the array
d_fraun = 2*(D^2)/lambda; %2*diagonal^2/lambda
d_NF = d_fraun/10; % the upper threshold to be in near-field region
d_bjo = 2*D; % the lower threshold to be in near-field region with constant amplitude
disp(['The near-field upper threshold is ' num2str(d_NF) ' (m)']);
disp(['Bjornson distance is ' num2str(d_bjo) ' (m)']);
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% The normal convention to evaluate position of the elements 
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end
az = pi/3;
el = pi/6;

%% Far-Field phase variation over surfce
x = UPA_Evaluate(lambda,M_V,M_H,az,el,d_V,d_H);
x = transpose(reshape(x,M_H,M_V));
x = flip(x,1);
figure;
imagesc(angle(x));
colormap(flipud(hot));
xlabel('$N_{\mathrm{H}}$','Interpreter','latex','FontSize',20);
ylabel('$N_{\mathrm{V}}$','Interpreter','latex','FontSize',20);
xticklabels([]);
yticklabels([]);
% Figure size
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);

%% Near-Field phase variation over surface
x = nearFieldChan(d_bjo,az,el,U,lambda);
x = transpose(reshape(x,M_H,M_V));
x = flip(x,1);
figure;
imagesc(angle(x));
colormap(flipud(hot));
xticklabels([]);
yticklabels([]);
xlabel('$N_{\mathrm{H}}$','Interpreter','latex','FontSize',20);
ylabel('$N_{\mathrm{V}}$','Interpreter','latex','FontSize',20);
% Figure size
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
