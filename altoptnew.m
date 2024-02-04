function [w,theta] = altoptnew(V,d)
N = length(d);

w = ones(N,1);
w = w/norm(w); % unit norm transmit beamformer

for i = 1:40
    theta = exp(1i*angle(w' * V))';
    w = (V*theta+d)/norm(V*theta+d);
end

end






