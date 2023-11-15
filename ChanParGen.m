function [azimuth,elevation,phase,dist] = ChanParGen(x,y,z,AP_coor,lambda)
% To generate the parameters for Saleh-Valenzuela channel model based on the 
% user location
%
%%%%%%%%%%%%%%%%%%  INPUT %%%%%%%%%%%%%%%%
%   x,y,z: are the coordinates of the UE in time
%   AP_coor : coordinate of the AP/BS/RIS 
%   lambda : wavelength
%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%
%   azimuth,elevation : Angle of arrivals at BS/RIS/AP
%   phase : phase shift due to the propagation delay
%   dist : distance to BS/RIS/AP 
 
d_x = x - AP_coor(1);
d_y = y - AP_coor(2);
d_z = z - AP_coor(3);
dist = sqrt(d_x.^2+d_y.^2+d_z.^2);
phase = exp(-1i*2*pi*dist/lambda);
azimuth = d_y./d_x;
azimuth = atan(azimuth);
elevation =  d_z ./ dist;
elevation = asin(elevation);

end