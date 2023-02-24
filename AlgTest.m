% This code is to test the functionality of the implemented algorithm to
% estimate the near field channel 
clc;clear;close all;
%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
%UPA Element configuration
M_H = 32; M_V = 32; M = M_H*M_V;
d_H = 1/2; d_V = 1/2; %In wavelengths
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D = sqrt(Hsize^2+Vsize^2); % Diagonal size of the array
d_fraun = 2*(D^2)/lambda/10; %2*diagonal^2/lambda/10
d_bjo = 2*D;
disp(['Fraunhofer distance is ' num2str(d_fraun) ' (m)']);
disp(['Bjornson distance is ' num2str(d_bjo) ' (m)']);
if d_fraun < d_bjo
    disp('The code will not work since bjornson distance is greater than fran dist');
    return;
end
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% The normal convention to evaluate position of the elements 
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end

%% Channel Estimation Parameters
% search resolution (It is very important)
varphiSRes = 2*M_H;
thetaSRes = 2*M_V;
distSRes = 4*M_H; % Distance resolution is very important to avoid convergence
Plim = M_H; % number of pilots

%Set the SNR
SNRdB_pilot = 10;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = 0;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H);
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 10;
nbrOfNoiseRealizations = 5;


%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(M-1,nbrOfAngleRealizations,nbrOfNoiseRealizations);
%% Iniitilize channel estimation 
% get the distance resolution for desinging the codebook
thre = 0.5; % correlation threshold
dmin = d_bjo;
dmax = 5*d_fraun;
dsearch = HeuristicDistRes(U,dmin,dmax,M,D,lambda,thre);

%Create a uniform grid of beams (like a DFT matrix) to be used at RIS
ElAngles = asin((-M_V/2:1:M_V/2-1)*2/M_V);
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
varphi_range = linspace(-pi/2,pi/2,varphiSRes);
theta_range = linspace(-pi/2,pi/2,thetaSRes);
dist_range = linspace(d_bjo,d_fraun,distSRes);
a_range = zeros(M,varphiSRes,thetaSRes,distSRes); % [M,Azimuth,Elevation,distance]
% obtain the array response vectors for all azimuth-elevation-distance
% triplet
for l =1:length(dist_range) % for each distance
    d_t = repelem(dist_range(l),1,varphiSRes);
    parfor i = 1:length(theta_range) % for each elevation
        a_range(:,:,i,l) = ...
            nearFieldChan(d_t,varphi_range,repelem(theta_range(i),1,varphiSRes),U,lambda);
    end
end

for n1 = 1:nbrOfAngleRealizations
    disp(n1);
    
    % generate the near-field channel 
    d_t = unifrnd(d_bjo,d_fraun,1);
    azimuth = unifrnd(-pi/3,pi/3,1);
    elevation = unifrnd(-pi/3,pi/3,1);
    g = nearFieldChan(d_t,azimuth,elevation,U,lambda); 

    var_amp_d= 64;
    d = sqrt(var_amp_d/2) * (randn + 1i*randn);

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);

    for n2 = 1:nbrOfNoiseRealizations
        disp(n2);
        % Estimate the channel using either all pilots or some of them
        %Select which two RIS configurations from the grid of beams to start with
        utilize = false(CBL,1);
        utilize(round(CBL/3)) = true;
        utilize(round(2*CBL/3)) = true;
        

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
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
     
            %Find an extra RIS configuration to use for pilot transmission
            if itr < Plim 

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

            end
        end
    end

end

rate = mean(mean(rate_proposed,3),2);

plot(2:Plim,repelem(mean(capacity),1,Plim-1));
hold on;
plot(2:Plim,rate(1:Plim-1));
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
fig.Children.TickLabelInterpreter = 'latex';
