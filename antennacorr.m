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
az = [0,pi/3];
corr = zeros(length(az),length(M_H));
for z = 1:length(az)
    azimuth = az(z);
    elevation = 0;
    d1 = 5;
    d2 = 15;
    % Evaluate the correlation between two distance for different array
    % apperture size 
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
        corr(z,l) = abs(a1'*a2/M(l))^2;
    end
end

% Plot the results
figure("defaultAxesFontSize",20,"DefaultLineLineWidth", 2,...
    "defaultAxesTickLabelInterpreter","latex");
for j = 1:length(az)
    plot(M,corr(j,:),'-o');
    hold on;
end
xlabel("Array Size M",'FontSize',20,'Interpreter','latex');
ylabel("Correlation",'FontSize',20,'Interpreter','latex');
xlim([0 max(M)]);
grid on;
legend1 = strcat("$\varphi = $",num2str(az(1)));
legend(legend1,"$\varphi = \pi /6$","Interpreter","latex");
title(("$d_1 =$ " + num2str(d1) + ", $d_2 =$ " + num2str(d2)),'Interpreter','latex');
