clear;clc;close all;
rng(7); % either 0 or 7
%% Room Scenario Initialization
Xmax = 5; % Room size = [-Xmax, Xmax]
Ymax = 5; % Room size = [-Ymax, Ymax]
numUE = 1; % number of users
RWL = 1000;  % Length of the random walk
Speedlim = [0.2,0.5]; % meter per second movement
RIS_center = [-Xmax,0,1.6]; % The RIS center point
z_t = repelem(1.5,numUE,RWL+1); % Z-coordinate of the user(does not change)

%% Rx/Tx/RIS Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 32; M_V = 32; M = M_H * M_V;
d_H_RIS = 1/2;d_V_RIS = 1/2;
Hsize = M_H*d_H_RIS*lambda;
Vsize = M_V*d_V_RIS*lambda;
disp(['Y-axis size is ' num2str(Hsize) ' (m)']);
disp(['Z-axis size is ' num2str(Vsize) ' (m)']);
d_fraun = 2*(Hsize^2+Vsize^2)/lambda/10; %2*diagonal^2/lambda/10
disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas
% This one is the normal convention to evaluate position of the elements
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H_RIS*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V_RIS*lambda; % dislocation with respect to center in z direction
    U(:,m) = [RIS_center(1); RIS_center(2)+ym; RIS_center(3)+zm]; %Position of the m-th element
end
%% 
d_BSRIS = 30;
% channel parameters (LOS mmWave) alpha + 10*beta*log10(d)
alpha = 61.4; % path loss reference in 1 (m) distance
beta = 1.46;
plEval = @(d) alpha + 10*beta*log10(d);
% Transmission params
noisepow = -96; % [dBm]
txpow = 0; % [dBm]

%% Plot configuration 
plt = true; % To plot the trajectory
pltconf = 'continous'; % 'continous' or 'discrete'

%% Random walk, get the farfield parameters
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Speedlim);
[azimuth,elevation,Cph,d_t] = ChanParGen(x_t,y_t,z_t,RIS_center,lambda);
if plt
    plotTrajectory(x_t,y_t,azimuth,elevation,pltconf,Xmax,Ymax,RIS_center);
    plotNearFieldBorder(x_t,y_t,Xmax,Ymax,RIS_center,d_fraun);
end
% report the near-field percentage
disp(['How much percentage in near field: ' num2str(sum(d_t < d_fraun,'all')/RWL/numUE*100)]);
% generate conventional and farfield channel
nearChan = nearFieldChanGen(x_t,y_t,z_t,U,lambda);
farChan = Cph .* UPA_Evaluate(lambda,M_V,M_H,azimuth,elevation,d_V_RIS,d_H_RIS);
% to understand how much we lose if we approximate the near field with far
% field
powgain = diag(abs(farChan'*nearChan).^2);
SE = log2(1+powgain);
figure('defaultLineLineWidth',2,'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',20);
plot(SE,'Color','b'); grid on; xlabel('Time (s)','Interpreter','latex');
ylabel('Spectral Efficiency (bit/s/Hz)','Interpreter','latex');
hold on; xlim([0,length(x_t)]);ylim([0,log2(1+max(powgain))+1]);
neart = d_t < d_fraun;
[~,idx] = find(neart==1);
plot(idx,SE(neart),'Marker','square','MarkerSize',10,'Color','r','LineStyle','none','LineWidth',2);
fig = gcf;
set(fig,'position',[60 50 1600 800]); % [left bottom width height]
