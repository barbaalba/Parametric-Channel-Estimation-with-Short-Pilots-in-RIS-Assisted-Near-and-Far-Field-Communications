function RIS_opt = altopt(H,g,d)
% Alternative Optimization function. It iteratively finds the optimum value
% for each RIS element
% Input: 
%   - H: channel between RIS and BS
%   - g: chanel betweeen RIS and UE
%   - d: direct channel between BS and UE
epoch = 20;
V = H*diag(g);  % cascaded channel
Vbar = V' * V;
N = length(g);
theta = ones(N,1);
rate = zeros(1,epoch);
% How many epoch
for j = 1:epoch
    % optimize over i-th RIS phase shift
    for i = 1:N
        theta(i) = exp(1i * angle(theta.'* Vbar(i,:).' + d.' * V(:,i)));
    end
    rate(j) = log2(1+(norm(V*theta + d)^2));
end
RIS_opt = theta;
plot(1:epoch,rate);
end