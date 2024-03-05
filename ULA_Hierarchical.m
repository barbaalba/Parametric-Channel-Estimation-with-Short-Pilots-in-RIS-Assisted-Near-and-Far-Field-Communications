function W = ULA_Hierarchical(N,d)
K = log2(N); % number of layer
W = zeros(N,2^K,K); %[Ant beams Layer]
for n = 1:2^K
    W(:,n,K) = ULA_Evaluate(N,asin(-1+(2*n-1)/N),d);
end
idx = zeros(2,1); 
for l = 1:K-1
    M = 2^(floor((l+1)/2)); % number of sub-arrays
    N_s = N/M; % number of antennas in sub-array
    N_A = M - M/2*mod(l,2); % number of active element per layer
    for m = 1:M
        idx(1) = (m-1)*N_s+1;
        idx(2) = m*N_s;
        if m <= N_A 
            W(idx(1):idx(2),1,K-l) = exp(-1i*m*(N_s-1)/N_s) * ...
                ULA_Evaluate(N_s,asin(-1+(2*m-1)/N_s),d);
        end
    end
    for n = 2:2^(K-l)
        W(:,n,K-l) = W(:,1,K-l) .* ULA_Evaluate(N,asin((2*(n-1))/2^(K-l)),d);
    end
end

end