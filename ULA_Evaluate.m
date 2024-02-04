function ULA_response = ULA_Evaluate(N,Azimuth,spacing)
% To evaluate normalized antenna vector for Uniform Linear Array 
% Inputs: 
%   lambda: Wavelength
%   N: Number of antenna array 
%   Azimuth: Azimuth of Arrival/Departure
% 
% Output:
%   ULA_response: Antenna response

d = spacing; % antenna spacing 

ULA_response = zeros(N,length(Azimuth));

for n = 1:length(Azimuth)

    ULA_response(:,n) = exp(-1i*2*pi*d*sin(Azimuth(n))*(0:N-1)');

end
