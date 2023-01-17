function phase = nearFieldChan(d_t,azimuth,elevation,U,lambda)
% To generate the real near field phase array without any approximation
% Input:
%   - d_t: the distance with respect to the reference element(center)
%   - U: is the relative location of each element from center of the array
% Output: 
%   - phase: is the phase vector when the full expression of near field
%            applied
U = -U;
M = size(U,2);
T = length(d_t);
d_m = zeros(M,T);
for t = 1:T
    d_m(:,t) = d_t(t)*sqrt(1+2*sin(azimuth(t))*cos(elevation(t))*U(2,:)/d_t(t)...
        + 2*sin(elevation(t))*U(3,:)/d_t(t) ...
        + (U(2,:).^2+U(3,:).^2) /d_t(t)^2 );
end
phase = exp(-1i*2*pi*d_m/lambda);
end