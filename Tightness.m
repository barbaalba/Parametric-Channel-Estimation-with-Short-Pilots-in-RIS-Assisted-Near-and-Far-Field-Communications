% This code is intended for validating the tightness of the 1st order
% approximation of teylor series for nearfield array response
clc;clear;close all;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 32; M_V = 32; M = M_H * M_V; % Size  of the arrat
d_H = 1/2;d_V = 1/2; % inter-element spacing
Ant_Center = [0,0,0]; 
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D = sqrt(Hsize^2+Vsize^2); % diagonal size of the array
d_fraun = 2*(D^2)/lambda/10; %2*diagonal^2/lambda/10
d_bjo = 2*D; % bjornson distance
disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
disp(['Bjornson distance is ' num2str(d_bjo) ' (m)']);
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas
% The relative position with respect to [0;0;0]
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end

T = 1000;
% Random realization of the channel
Azimuth = unifrnd(-pi/3,pi/3,1,T);
Elevation = unifrnd(-pi/3,pi/3,1,T);
dist = unifrnd(d_bjo,d_fraun,1,T);
nearcorr = zeros(1,T);
cap = zeros(1,T);
farcorr = zeros(1,T);
% Evaluate the correlation between the real channel and near-field and
% far-field approximation
parfor t = 1:T
    nearchan = nearFieldChan(dist(t),Azimuth(t),Elevation(t),U,lambda);
    cap(t) = abs(1/M*nearchan'*nearchan)^2;
    nearchanapp = nearFieldApproxChan(dist(t),Azimuth(t),Elevation(t),U,lambda);
    nearcorr(t) = abs(1/M*nearchan'*nearchanapp)^2;
    farchan = UPA_Evaluate(lambda,M_V,M_H,Azimuth(t),Elevation(t),d_V,d_H);
    farcorr(t) = abs(1/M*farchan'*nearchan)^2;
end

% plot
figure("defaultAxesFontSize",20,"DefaultLineLineWidth", 2,...
    "defaultAxesTickLabelInterpreter","latex");
plot(pow2db(cap),'k');
hold on;
plot(pow2db(nearcorr));
plot(pow2db(farcorr));
xlabel("Channel Instances","FontSize",20,Interpreter="latex");
ylabel("Channel Correlation","FontSize",20,"Interpreter","latex");
legend("Real Channel","1st order Approx","Far-Field Approx",'interpreter','latex')