function [capacity,R,idx] = MLEfunction3D(G,h,hd,U,dsearch,M_H,M_V,d_H,d_V,lambda,SNR_pilot,SNR_data,Prep)
% This function estimates the phase and angles and distance based on MLE
dsearch = flip(dsearch); % to go from far-field to near-field distances
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
ElAngles = asin(-1:2/M_V*length(dsearch):1);
AzAngles = asin((-M_H/2:1:M_H/2-1)*2/M_H);
CBL = length(ElAngles)*length(AzAngles)*length(dsearch); % Code book length
beamPar = zeros(CBL,3); % Elevation-Azimuth pair
beamresponses = zeros(M,CBL); %[Ant,Code book size]
% The code book parameters. Each row is one beamforming vector [Az,El,dist]
for l=1:length(dsearch)
    d_t = repelem(dsearch(l),1,length(AzAngles));
    for i = 1:length(ElAngles)
        lowrange = (i-1)*M_H+(l-1)*length(ElAngles)*M_H+1; % indexing variable
        uprange = i*M_H+(l-1)*length(ElAngles)*M_H; % indexing variable
        beamPar(lowrange:uprange,:,:) = ...
            [AzAngles.' , repelem(ElAngles(i),1,length(AzAngles)).' , d_t.'];
    end
end
% Build the dictionary array response
parfor i = 1:CBL
    beamresponses(:,i) = nearFieldChan(beamPar(i,3), beamPar(i,1), beamPar(i,2),U,lambda);
end

% Define a fine grid of angle directions and distance to be used in
% estimator
varphi_range = linspace(-pi/2,pi/2,SRes);
theta_range = linspace(-pi/2,pi/2,SRes);
dist_range = linspace(0.1,5.1,SRes);
a_range = zeros(M,SRes,SRes,SRes); % [M,Azimuth,Elevation,distance]
% obtain the array response vectors for all azimuth-elevation-distance
% triplet
for l =1:SRes % for each distance
    d_t = repelem(dist_range(l),1,SRes);
    parfor i = 1:SRes % for each elevation
        a_range(:,:,i,l) = ...
            nearFieldChan(d_t,varphi_range,repelem(theta_range(i),1,SRes),U,lambda);
    end
end

idx = zeros(nbrOfAngleRealizations,nbrOfNoiseRealizations);
n1 = 1;
while n1<=nbrOfAngleRealizations
    disp(n1);
    
    for n2 = 1:nbrOfNoiseRealizations
        % Estimate the channel using either all pilots or some of them
        %if n1 == 1
            %Select which two RIS configurations from the grid of beams to start with
            utilize = false(CBL,1);
            if n1 == 1
                % For the first time we initialize randomly 
                utilize(round(CBL/3)) = true;
                utilize(round(2*CBL/3)) = true;
                Plim = M; % number of pilots
            else
                % Start the estimation with using the previous best RIS
                % config to track the channel
                utilize(bestInit1) = true;
                idx(n1,n2) = randi(CBL-M_H)+M_H;
                while idx(n1,n2) == bestInit1
                    idx(n1,n2) = randi(CBL-M_H)+M_H;
                end
                utilize(idx(n1,n2)) = true;
                Plim = M; % number of pilots
            end
            g = G(:,n1);
            d = hd(n1);
            %Compute the exact capacity for the system Eq. (3)
            capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);
            %Define the initial transmission setup
%             RIS_directions = beamPar(utilize,:); % [Az,El,dist]
            RISconfigs = Dh_angles*beamresponses(:,utilize);
            B = RISconfigs';
            %Generate the noise
            noise = (randn(CBL,1)+1i*randn(CBL,1))/sqrt(2);

            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1

                %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
                
                % Estimate the Channel using the developed MLE
                [~,var_phas_d_est,~,~,g_est,~,~] = MLE3D(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);
             

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
                    %RIS_directions = beamAngles(utilize,:);
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