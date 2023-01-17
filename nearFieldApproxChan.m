function phase = nearFieldApproxChan(d_t,azimuth,elevation,U,lambda)
% This function use the first order approximation of taylor series to
% approximate the near field 
% Input
%   - d_t: is the distance for the reference element
%   - U: is the relative distance of each element from the array's center
% Output
%   - phase: The phase of all array considering the near field approx
U = -U;
M = size(U,2);
T = length(d_t);
d_m = zeros(M,T);
for t = 1:T
   d_m(:,t) = d_t(t) + sin(azimuth(t))*cos(elevation(t))*U(2,:) + ...
       sin(elevation(t))*U(3,:) + U(2,:).^2/(2*d_t(t)) + U(3,:).^2/(2*d_t(t));
end

phase = exp(-1i*2*pi*d_m/lambda);

end