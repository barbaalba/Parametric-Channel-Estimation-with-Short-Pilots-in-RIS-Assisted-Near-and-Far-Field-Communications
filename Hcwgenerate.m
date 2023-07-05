function W = Hcwgenerate(N,d,lambda)
% this function is constructed following Z. Xiao, T. He, P. Xia and 
% X. -G. Xia, "Hierarchical Codebook Design for Beamforming Training in 
% Millimeter-Wave Communication," in IEEE Transactions on Wireless 
% Communications, vol. 15, no. 5, pp. 3380-3392, May 2016

K = log2(N); % number of layers in hierarchical design
W = struct;
for i = 0:1:K
    eval(['W.Layer' num2str(i) '= [];']);
end
fl = fieldnames(W);
% The last layer which has sharpest beam
n = 1:N; azimuth = asin(-1+(2*n-1)/N);
W.(fl{end}) = ULA_Evaluate(lambda,N,azimuth,d);

for l = 1:K
    k = K-l;
    M = 2^floor((l+1)/2); % subarrays
    N_s = N/M; % number of subarray elements
    N_a = (1-mod(l,2))*M + mod(l,2)*M/2; % subarray activation
    m = 1:N_a; azimuth = asin(-1+(2*m-1)/N_s);
    f_m = exp(-1i*m*(N_s-1)/N_s*pi).*ULA_Evaluate(lambda,N_s,azimuth,d);
    % The first beam in each layer
    w = zeros(N,2^k); 
    w(1:N_a*N_s,1) = reshape(f_m,N_s*N_a,1);
    n = 2:2^k; azimuth = asin(2*(n-1)/N);
    w(:,2:end) = w(:,1) .* ULA_Evaluate(lambda,N,azimuth,d);
    W.(fl{k+1}) = w;
end

end