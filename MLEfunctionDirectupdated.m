function [capacity,R,idx] = MLEfunctionDirectupdated(G,h,hd,M_H,M_V,d_H,d_V,lambda,SNR_pilot,SNR_data,Prep)
% This function estimates the phase and angles and then design the best RIS
% configuration. The algorithm is based on: M. Haghshenas, P. Ramezani
% E. Bjornson, "Efficient LOS Channel Estimation for RIS-Aided 
% Communications Under Non-Stationary Mobility", ICC 2023

Dh = diag(h);
Dh_angles = diag(h./abs(h)); 
M = M_H*M_V;
nbrOfAngleRealizations = size(G,2);
nbrOfNoiseRealizations = 1;
SRes = 2*M_H; % search resolution

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(nbrOfAngleRealizations,nbrOfNoiseRealizations);

%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
ElAngles = asin((-M_V/2:1:M_V/2-1)*2/M_V);
AzAngles = asin((-M_H/2:1:M_H/2-1)*2/M_H);
beamAngles = zeros(M,2); % [El,Az]
% Row i in beamAngles corresponds to the column i in beamresponsess
beamresponses = zeros(M,M);
for i = 1:length(ElAngles)
    beamAngles ((i-1)*length(AzAngles)+1: i*length(AzAngles),:) = [repelem(ElAngles(i),length(AzAngles),1) AzAngles'];
    beamresponses(:,(i-1)*length(AzAngles)+1: i*length(AzAngles)) = UPA_Evaluate(lambda,M_V,M_H,AzAngles,repelem(ElAngles(i),1,length(AzAngles)),d_V,d_H);
end

% Define a fine grid of angle directions to analyze when searching for angle of arrival
varphi_range = linspace(-pi/2,pi/2,SRes);
theta_range = linspace(-pi/2,pi/2,SRes);
a_range = zeros(M,SRes,SRes); % [M,Azimuth,Elevation]
% obtain the array response vectors for all azimuth-elevation pairs
parfor i = 1:SRes
    a_range(:,:,i) = ...
    UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,SRes),d_V,d_H);
end

idx = zeros(nbrOfAngleRealizations,nbrOfNoiseRealizations);
n1 = 1;
while n1<=nbrOfAngleRealizations
    disp(n1);
    for n2 = 1:nbrOfNoiseRealizations
        % Estimate the channel starting with n1 symbol
        %if mod(n1,Prep) == 1
            %Select which two RIS configurations from the grid of beams to start with
            utilize = false(M,1);
            if n1 == 1
                % For the first time we initialize randomly 
                utilize(round(M/3)) = true;
                utilize(round(2*M/3)) = true;
                Plim = M; % number of pilots
            else
                % Start the estimation with using the previous best RIS
                % config to track the channel
                utilize(bestInit1) = true;
                idx(n1,n2) = randi(M-M_H)+M_H;
                while idx(n1,n2) == bestInit1
                    idx(n1,n2) = randi(M-M_H)+M_H;
                end
                utilize(idx(n1,n2)) = true;
                Plim = M; % number of pilots
            end
            g = G(:,n1);
            d = hd(n1);
            %Compute the exact capacity for the system Eq. (3)
            capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);
            %Define the initial transmission setup
            RISconfigs = Dh_angles*beamresponses(:,utilize);
            B = RISconfigs';
            %Generate the noise
            noise = (randn(M,1)+1i*randn(M,1))/sqrt(2);

            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1
                %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
                
                % Estimate the Channel using the developed MLE
                [~,var_phas_d_est,~,~,g_est,~,~] = MLE(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);
             
             
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig = angle(Dh*g_est)-var_phas_d_est;
     

                %Compute the corresponding achievable rate
                rate_proposed(n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);


                %Find an extra RIS configuration to use for pilot transmission
                if itr < Plim-1

                    %Find which angles in the grid-of-beams haven't been used
                    unusedBeamresponses = beamresponses;
                    unusedBeamresponses(:,utilize==1) = 0;

                    %Guess what the channel would be with the different beams
                    guessBeam = Dh*unusedBeamresponses;
               


                    %Find which of the guessed channels matches best with the
                    %currently best RIS configuration

                    closestBeam = abs(exp(-1i*RISconfig).'*guessBeam); 
                    [~,bestUnusedBeamidx] = max(closestBeam);


                    %Add a pilot transmission using the new RIS configuration
                    utilize(bestUnusedBeamidx) = true;
                    RISconfigs = Dh_angles*beamresponses(:,utilize);
                    B = RISconfigs';
                else
                    closestBeam = abs(exp(-1i*RISconfig).'*Dh_angles*beamresponses); 
                    [~,bestInit1] = max(closestBeam);
                    % To correct the drops due to the pilots during channel
                    % esitimation
                    rate_proposed(n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
                    n1 = n1+1;
                end
            end
        %else
%             g = G(:,n1);
%             d = hd(n1);
%             capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d))^2);
%             rate_proposed(n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
%             n1 = n1+1;
        %end
    end
end
R = mean(rate_proposed,2);
end