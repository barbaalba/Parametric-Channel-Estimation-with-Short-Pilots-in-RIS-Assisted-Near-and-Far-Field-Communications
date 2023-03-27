% This code is written to compare the far-field and near-field eignvalues 
% to confirm that both near-field and far-field lies on the subspace
% with equal dimension
clear;clc;close all;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
% Antenna configuration
d_H = 1/2;d_V = 1/2;
M_H = 32; M_V = M_H; M = M_H*M_V;
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

corrNear = zeros(M,M);
corrFar = zeros(M,M);
parfor i = 1:M_H
    az = unifrnd(-pi/2,pi/2);
    el = unifrnd(-pi/2,pi/2);
    HNear = nearFieldChan(d_bjo,az,el,U,lambda);
    Hfar = UPA_Evaluate(lambda,M_V,M_H,az,el,d_V,d_H);
    corrNear = corrNear + HNear*HNear';
    corrFar = corrFar + Hfar*Hfar';
end
corrNear = corrNear/M;
corrFar = corrFar/M;
eNear = flip(eig(corrNear));
eNear = eNear(eNear>0);
eFar = flip(eig(corrFar));
eFar = eFar(eFar>0);
plot(pow2db(eNear),'LineWidth',4);
hold on
plot(pow2db(eFar),'LineWidth',2,'LineStyle','--');
xlabel('Ordered Eigenvalues','FontSize',20,'Interpreter','latex');
ylabel('Eigenvalues [dB]','FontSize',20,'Interpreter','latex');
legend('Near Field','FarField','FontSize',20,'interpreter','latex');
ax = gca;
ax.FontSize = 20;
grid on;
title("RIS with size of " +M_H +"$\times$" +M_V+" and element spacing of "+d_H,'Interpreter','latex')