function [d_est,g_est] = MLE4D(y,L,F,a_range,SNR_pilot)
M = length(y)/L; % BS antenna

thetaSRes = size(a_range,3);
varphiSRes = size(a_range,2);
distSRes = size(a_range,4);
%Compute the ML utility function for all potential AOA/dist triplet
utility_num = zeros(varphiSRes,thetaSRes,distSRes); %numerator of the utility function
utility_den = zeros(varphiSRes,thetaSRes,distSRes); %denominator of the utility function
% [Az,El,dist]
for l = 1:distSRes
    for i = 1:thetaSRes
        utility_num(:,i,l) = abs( y' * ( eye(M*L) - ( kron(ones(L,L),eye(M)) ) /L ) ...
            * F * a_range(:,:,i,l) ).^2;

        utility_den(:,i,l) = sum(abs(F * a_range(:,:,i,l)).^2,1) - 1/L * ...
            sum(abs(kron(ones(L,1).',eye(M)) * F * a_range(:,:,i,l)).^2,1);
    end
end

utilityfunction = utility_num ./ utility_den; %[Az,El,dist]
%Extract the angle estimate
[~,maxind] = max(utilityfunction,[],'all');
[Azidx,Elidx,didx] = ind2sub([varphiSRes,thetaSRes,distSRes],maxind);
a = a_range(:,Azidx,Elidx,didx);

%Estimate g
var_amp_g_num = abs(y' * ( eye(M*L) - ( kron(ones(L,L),eye(M)) ) /L ) * F * a);
var_amp_g_den = sum(abs(F * a).^2,1) - 1/L * ...
    sum(abs(kron(ones(L,1).',eye(M)) * F * a).^2,1);
var_amp_g = 1/sqrt(SNR_pilot) * var_amp_g_num / var_amp_g_den;
phase_g = - angle(y' * ( eye(M*L) - ( kron(ones(L,L),eye(M)) ) /L ) * F * a);
g_est = var_amp_g * exp(1i*phase_g) * a;
d_est = 1/(sqrt(SNR_pilot) * L) * kron(ones(L,1).',eye(M)) * (y - sqrt(SNR_pilot)*F*g_est);
end

