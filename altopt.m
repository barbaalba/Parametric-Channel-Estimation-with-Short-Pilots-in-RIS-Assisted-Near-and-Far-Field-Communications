function RIS_opt = altopt(H,g,d)
% Alternative Optimization function. It iteratively finds the optimum value
% for each RIS element
% Input: 
%   - H: channel between RIS and BS
%   - g: chanel betweeen RIS and UE
%   - d: direct channel between BS and UE
epoch = 20;
G = H*diag(g);  % cascaded channel
Gbar = G' * G;
N = length(g);
theta = ones(N,1);
rate = zeros(1,epoch);

% How many epoch
for j = 1:epoch
    % optimize over i-th RIS phase shift
    for i = 1:N
        idx = 1:N;
        idx(i) = [];
        theta(i) = exp(-1i * angle(theta(idx)'* Gbar(idx,i) + d' * G(:,i)));
    end
    rate(j) = log2(1+(norm(G*theta + d)^2));
end
RIS_opt = theta;
plot(1:epoch,rate);
end