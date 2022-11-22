function [UPA_respones,U] = UPA_Evaluate(lambda,M_V,M_H,Azimuth,Elevation,vd,hd)
%   Evaluate the planar antenna responses
%
%   Inputs:
%       lambda: Wavelength
%       M_V:    Number of vertical antenna elements
%       M_H:    Number of horizental antenna elements
%       Azimuth
%       Elevation
%
%   Outputs:
%       UPA_response: Antenna response with respect to Azimuth and
%       Elevation
%       U:  Antenna Elements positions
%           Arrays are on Y-Z Axis. The index of the array is as follows:
%           6 7 8 
%           3 4 5
%           0 1 2
%

%% Parameters Initializations

d_H = hd*lambda; % Horizontal antenna spacing    
d_V = vd*lambda; % Vertical antenna spacing

M = M_H*M_V; % Total number of antennas
U = zeros(3,M); % Matrix containing the position of the antennas

i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index

%% Element position evaluation

% This one is the normal convention to evaluate position of the elements
for m = 1:M
    
        U(:,m) = [0; i(m)*d_H; j(m)*d_V]; %Position of the m-th element

end

%% Array Response Evaluation
%{ 
The evaluation is based on the following book:
Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency
by Emil Björnson, Jakob Hoydis, Luca Sanguinetti
%}

UPA_respones = zeros(M,length(Azimuth));

for n = 1:length(Azimuth)
    
    UPA_respones(:,n) = functionSpatialSignature3DLoS(U,Azimuth(n),...
        Elevation(n),lambda);

end
