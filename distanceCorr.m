% This code is to understand better the disance correlation, we notice that
% by getting further from the antenna array, the correlation starts to
% increase meaning that angles are the major parameters to determine the
% beam
clear; clc; close all;
%% Initilization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 64; M_V = 64; M = M_H * M_V;
d_H_RIS = 1/2;d_V_RIS = 1/2;
Hsize = M_H*d_H_RIS*lambda;
Vsize = M_V*d_V_RIS*lambda;
d_fraun = 2*(Hsize^2+Vsize^2)/lambda/10; %2*diagonal^2/lambda/10
disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
%% Evaluate the position of the elements with respect to the center [0;0;0] 
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% Near Field approximation based on distance and direction of arrivals
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H_RIS*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V_RIS*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end
%% Average the correlation over different angles at specific distances
T = 1000;
azimuth = rand(1,T)*pi-pi/2;
elevation = rand(1,T)*pi-pi/2;
d1 = repelem(0.1,1,T);
d2 = d1:0.1:d_fraun;
dist = d2 - d1(1);
d2 = d2';
d2 = repmat(d2,1,T);
corr = zeros(1,size(d2,1));
a1 = nearFieldChan(d1,azimuth,elevation,U,lambda); 
parfor i = 1:size(d2,1)
a2 = nearFieldChan(d2(i,:),azimuth,elevation,U,lambda);
corr(i) = mean(abs(diag(a1'*a2)/M).^2);
end
% plot the correlation
plot(dist,pow2db(corr));