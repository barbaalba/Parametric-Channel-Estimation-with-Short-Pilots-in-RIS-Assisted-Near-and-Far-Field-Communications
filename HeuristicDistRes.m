function dsearch = HeuristicDistRes(U,dmin,dmax,M,D,lambda,thre)
% this function produce the distance resolutions Heuristically for MLE.
% We use these distances to design code book 
% Input: 
%   - U: the relative poistion of the RIS element with respect to [0;0;0]
%   - dmin: the minimum distance that user might be 
%   - dmax: the maximum distance with respect to user location
%   - M: number of total elements
%   - D: is the diagonal length of array apperture
%   - thre: is the correlation threshold with which we decide on the next
%   distance resolution
% Output:
%   - The distances that we use to design codebook 
d1 = dmin; d2list=d1:D/10:dmax;
l = 1;
dsearch = zeros(1,length(d2list)); % initialize the final distance search
dsearch(1) = d1;
while l < length(d2list)+1
    d2 =d2list(l);
    dist = 0.5/d1 - 0.5/d2;
    corr = 1/(M^2) * abs(sum(exp(2i*pi/lambda*dist*(U(2,:).^2+U(3,:).^2))))^2;
    if corr < thre
        d1 = d2;
        dsearch(l) = d2;
    end
    l = l+1;
end
dsearch = dsearch(dsearch>0);

end