function W = UPA_Hierarchical(N_H,N_V,d_H,d_V)
W_H = ULA_Hierarchical(N_H,d_H);
W_V = ULA_Hierarchical(N_V,d_V);
N = N_H * N_V;
W = zeros(N,N,size(W_H,3));
for k = 1:size(W_H,3)
    beam = zeros(N,4^k);
    beam_H = W_H(:,1:2^k,k);
    beam_V = W_V(:,1:2^k,k);
    for i = 1:2^k
        beam(:,(i-1)*2^k+1:2^k*i) = kron(beam_V,beam_H(:,i));
    end
    W(:,1:4^k,k) = beam;
end
end