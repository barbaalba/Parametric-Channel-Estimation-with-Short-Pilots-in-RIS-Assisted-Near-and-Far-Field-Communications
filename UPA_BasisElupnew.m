function [theta,phi,bsnum] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,azref,elref)
% This function evaluate the orthogonal basis for UPA structure depending
% on the initial points
% Input 
%   M_V : number of antennas in vertical axis
%   M_H : number of antennas in horizental axis
%   d_V : normalized antenna spacing in vertical axis
%   d_H : normalized antenna spacing in horizental axis
%   azref : the reference azimuth angle
% Output
%   theta,phi : angle pairs
%   bsnum : number of orthogonal beams

L_V = M_V * d_V;
L_H = M_H * d_H;

%% first find the elevation orthogonal to 90 degree
theta = []; % 0 as the elevation reference

for k = -2*d_V*M_V:1:2*d_V*M_V
    omega = k / L_V;
    val = asin(omega + sin(elref));
    if isreal(val)
        if (val == pi/2 || val == -pi/2) && d_V == 1/2 && ismember(-val,theta)
        
        else
            theta = [theta val];
        end
    end
end

%% Initialize the azimuth angles with the reference angle
phi = struct;
for j = 1:length(theta)
    eval(['phi.theta' num2str(j) '= [];']);
end
fn = fieldnames(phi);

%% Build the azimuth for each elevation angle
for j = 1:length(theta)
    for k = -2*M_H*d_H:1:2*d_H*M_H % -2::2
        gama = k/L_H;
        gama = gama/cos(theta(j)) + sin(azref);
        val = asin(gama);
        if isreal(val)
            if (val == pi/2 || val == -pi/2) && d_H == 1/2 && ismember(-val,phi.(fn{j}))

            else
                phi.(fn{j}) = [phi.(fn{j}) val];
            end
        end
    end
end
fn = fieldnames(phi);
bsnum = 0; % Basis number
for j = 1:length(fn)
    bsnum = bsnum + length(phi.(fn{j}));
end
disp(['Number of basis before process: ' num2str(bsnum)]);

end