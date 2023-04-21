function W = UPA_Codebook(lambda,theta,phi,M_V,M_H,d_V,d_H)
% Generate the orthogonal beams for an UPA 
% Input
%   - theta: Elevations angles that are orthogonal to each other
%   - phi: Azimuth angles that are orthogonal on the same elevation angles

fn = fieldnames(phi);
bsnum = 0; % Basis number
for j = 1:length(fn)
    bsnum = bsnum + length(phi.(fn{j}));
end
disp(['Number of basis after process: ' num2str(bsnum)]);
% Initialize the basis precoder
W = zeros(M_V*M_H,bsnum);
dumbval=1; % dummy variable
for j = 1:length(theta)
    for k = 1:length(phi.(fn{j}))
        W(:,dumbval) = UPA_Evaluate(lambda,M_V,M_H,phi.(fn{j})(k),theta(j),d_V,d_H);
        dumbval = dumbval+1;
    end
end
end