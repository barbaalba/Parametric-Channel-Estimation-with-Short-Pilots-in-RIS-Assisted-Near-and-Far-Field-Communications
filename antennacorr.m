% This code is intended to demonstrate how nearfield behave different than
% Far-field 
clear; clc; close all;
%%
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
% Antenna configuration
d_H = 1/2;d_V = 1/2;
M_H = 20:5:300;M_V = M_H; M = M_H.*M_V;
% Point setup in space
azimuth = pi/3;
elevation = 0;
d1 = 5;
d2 = 15;
% Evaluate the correlation between two distance for different array
% apperture size 
corr = zeros(1,length(M_H));
for l = 1:length(M_H)
    Hsize = M_H(l)*d_H*lambda;
    Vsize = M_V(l)*d_V*lambda;
    d_fraun = 2*(Hsize^2+Vsize^2)/lambda/10; %2*diagonal^2/lambda/10
    disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
    i = @(m) mod(m-1,M_H(l)); % Horizontal index
    j = @(m) floor((m-1)/M_H(l)); % Vertical index
    U = zeros(3,M(l)); % Matrix containing the position of the RIS antennas 
    % Near Field approximation based on distance and direction of arrivals
    for m = 1:M(l)  
        ym = (-(M_H(l)-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
        zm = (-(M_V(l)-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
        U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
    end
    a1 = nearFieldChan(d1,azimuth,elevation,U,lambda); 
    a2 = nearFieldChan(d2,azimuth,elevation,U,lambda);
    corr(l) = abs(a1'*a2/M(l))^2;
end
plot(M,corr,'-o');