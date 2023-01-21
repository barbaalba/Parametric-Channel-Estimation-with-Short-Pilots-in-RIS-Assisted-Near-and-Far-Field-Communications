clear;clc;close all;
%rng(9); % either 9 or 7
%% Room Scenario Initialization
Xmax = 1.5; % Room size = [-Xmax, Xmax]
Ymax = 1.5; % Room size = [-Ymax, Ymax]
Xinit = -1.3;Yinit = 0; % Initial location of the user in the room
numUE = 1; % number of users
RWL = 20;  % Length of the random walk
Speedlim = [0.1,0.4]; % meter per second movement
RIS_center = [-Xmax,0,1.6]; % The RIS center point
z_t = repelem(1.5,numUE,RWL+1); % Z-coordinate of the user(does not change)

%% Rx/Tx/RIS Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

% Set the SNR
SNRdB_pilot = 0;
SNR_pilot = db2pow(SNRdB_pilot);
SNRdB_data = SNRdB_pilot-10;
SNR_data = db2pow(SNRdB_data);
Prep = 1;

% RIS parameters
M_H = 16; M_V = 16; M = M_H * M_V;
d_H_RIS = 1/2;d_V_RIS = 1/2;
Hsize = M_H*d_H_RIS*lambda;
Vsize = M_V*d_V_RIS*lambda;
disp(['Y-axis size is ' num2str(Hsize) ' (m)']);
disp(['Z-axis size is ' num2str(Vsize) ' (m)']);
D = sqrt(Hsize^2+Vsize^2);
d_fraun = 2*(D^2)/lambda/10; %2*diagonal^2/lambda/10
d_bjo = 2*D;
disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
disp(['Bjornson distance is ' num2str(d_bjo) ' (m)']);
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% The normal convention to evaluate position of the elements in the room
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H_RIS*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V_RIS*lambda; % dislocation with respect to center in z direction
    U(:,m) = [RIS_center(1); RIS_center(2)+ym; RIS_center(3)+zm]; % Position of the m-th element
end
%% Known channel between RIS and RSU/BS
varphi_BS = pi/6;
theta_BS = 0;
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V_RIS,d_H_RIS);

%% The channel between UE and BS
var_amp_d= 64;
hd = sqrt(var_amp_d/2) * (randn(1,RWL+1) + 1i*randn(1,RWL+1)); 

%% Plot configuration 
plt = false; % To plot the trajectory
pltconf = 'continous'; % 'continous' or 'discrete'

%% Random walk, get the farfield parameters
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Xinit,Yinit,Speedlim);
[azimuth,elevation,Cph,d_t] = ChanParGen(x_t,y_t,z_t,RIS_center,lambda);
if plt
    plotTrajectory(x_t,y_t,azimuth,elevation,pltconf,Xmax,Ymax,RIS_center);
    plotNearFieldBorder(x_t,y_t,Xmax,Ymax,RIS_center,d_fraun);
end
% report the near-field percentage
disp(['How much percentage in near field: ' num2str(sum(d_t < d_fraun,'all')/RWL/numUE*100)]);
% generate the true and farfield channel
G = realChanGen(x_t,y_t,z_t,U,lambda); % it is the real channel
farChan = Cph .* UPA_Evaluate(lambda,M_V,M_H,azimuth,elevation,d_V_RIS,d_H_RIS);
% The relative position with respect to [0;0;0]
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H_RIS*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V_RIS*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end

% Near field approximation based on distance and angles
nearChanapprox = nearFieldApproxChan(d_t,azimuth,elevation,U,lambda);
% Near field full expression
nearChan = nearFieldChan(d_t,azimuth,elevation,U,lambda); 

%% Real Channel V.s. Near Field V.s. Far Field approximation
% Maximum gain
powgain = diag(abs(G'*G).^2);
SE = log2(1+powgain);
figure('defaultLineLineWidth',2,'defaultAxesTickLabelInterpreter','latex','defaultAxesFontSize',20);
plot(SE,'Color','k','LineWidth',4); grid on; xlabel('Time (s)','Interpreter','latex');
ylabel('Spectral Efficiency (bit/s/Hz)','Interpreter','latex');
hold on; xlim([0,length(x_t)]);ylim([0,log2(1+max(powgain))+1]);
% Near Field approximation
powgain = diag(abs(nearChanapprox'*G).^2);
SE = log2(1+powgain);
plot(SE,'r');
% Near Field without approximation
powgain = diag(abs(nearChan'*G).^2);
SE = log2(1+powgain);
plot(SE,'g');
% Far Field approximation
powgain = diag(abs(farChan'*G).^2);
SE = log2(1+powgain);
plot(SE,'Color','b','LineStyle',':'); grid on; xlabel('Time (s)','Interpreter','latex');
fig = gcf;
set(fig,'position',[60 50 1600 800]); % [left bottom width height]
legend('Max Gain','Near Field Approx','Near Field','Far Field Approx','interpreter','latex');
title([num2str(M_H) '\times' num2str(M_V)]);
% Mark the detected near filed time instances
neart = d_t < d_fraun;
[~,idx] = find(neart==1);
%plot(idx,SE(neart),'Marker','square','MarkerSize',10,'Color','r','LineStyle','none','LineWidth',2);
%% Channel Estimation 

% get the distance resolution for desinging the codebook
thre = 0.7; % correlation threshold
dmin = min([d_bjo,d_t]);
dmax = max(d_t);
dsearch = HeuristicDistRes(U,dmin,dmax,M,D,lambda,thre);

[capacity,R,~] = MLEfunction3D(G,h,hd,U,dsearch,M_H,M_V,d_H_RIS,d_V_RIS,lambda,SNR_pilot,SNR_data,Prep);
