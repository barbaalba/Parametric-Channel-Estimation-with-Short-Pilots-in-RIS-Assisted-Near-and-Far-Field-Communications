function R = UPAcorrelation(N_H,N_V,d_H,d_V,lambda)
N = N_H*N_V;
i = @(m) mod(m-1,N_H); % Horizontal index
j = @(m) floor((m-1)/N_H); % Vertical index
U = zeros(3,N); % Matrix containing the position of the RIS antennas 
Dist = zeros(N,N);
% The normal convention to evaluate position of the elements 
for m = 1:N 
    ym =  i(m)*d_H*lambda; % dislocation with respect to center in y direction
    zm = j(m)*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end

Dist = pdist(U.');
Dist = squareform(Dist);
R = sinc(2/lambda*Dist); % the correlation fucntion between two elements
end