function [d_est,var_phas_d_est,var_amp_g_est,var_phas_g_est,g_est,Azidx,Elidx] = MLE(y,L,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot)
% Input:
%   - y: rx signal
%   - L: number of pilots
%   - 
SRes = size(a_range,3);
%Compute the ML utility function for all potential angle pairs
utility_num = zeros(SRes,SRes); %numerator of the utility function
utility_den = zeros(SRes,SRes); %denominator of the utility function
% [Az,El] Search
for i = 1:SRes
    utility_num(:,i) = abs(y' * (eye(L) - (L)^-1 * ones(L,L))* ...
        B*Dh*a_range(:,:,i)).^2;
    utility_den(:,i) = sum(abs(B*Dh*a_range(:,:,i)).^2,1) - (L)^-1 * ...
        abs(ones(1,L)*B*Dh*a_range(:,:,i)).^2;
end
           
             
utilityfunction = utility_num ./ utility_den;

%Extract the angle estimate
[~,maxind] = max(utilityfunction,[],'all');
[Azidx,Elidx] = ind2sub([SRes,SRes],maxind);
a = a_range(:,Azidx,Elidx);

%Estimate g
var_amp_g_num = abs(y' * (eye(L) - (L)^-1 * ones(L,L))* B*Dh*a)^2;
var_amp_g_den = SNR_pilot * (sum(abs(B*Dh*a).^2,1) - ...
    (L)^-1 * abs(ones(1,L)*B*Dh*a)^2)^2;
var_amp_g_est = var_amp_g_num/var_amp_g_den;
var_phas_g_est = - angle (y' * (eye(L) - (L)^-1 * ones(L,L))* B*Dh*a);
g_est = sqrt(var_amp_g_est) * exp(1i*var_phas_g_est) *a;

%Estimate d
var_amp_d_est = (SNR_pilot)^-1 * (L)^-2 * abs (ones(1,L) * (y - sqrt(SNR_pilot)*B*Dh*g_est))^2;
var_phas_d_est = angle (ones(1,L) * (y - sqrt(SNR_pilot)*B*Dh*g_est));
d_est = sqrt(var_amp_d_est) * exp(1i*var_phas_d_est);

end