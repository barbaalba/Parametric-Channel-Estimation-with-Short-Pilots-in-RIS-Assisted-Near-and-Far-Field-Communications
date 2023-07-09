% This code is to test the functionality of the implemented algorithm to
% estimate the channel by using widebeams
clc;clear;close all;
%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

%UPA Element configuration
M_H = 16; M_V = 16; M = M_H*M_V;
d_H = 1/2; d_V = 1/2; %In wavelengths

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
varphiSRes = 8*M_H;
thetaSRes = 8*M_V;
Plim = 32; % number of pilots

%Set the SNR
SNRdB_pilot = 0;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -10;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H);
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 100;
nbrOfNoiseRealizations = 5;

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_wide = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
SNR_total = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
SNR_w_total = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);

%% Iniitilize channel estimation 
% Generate code book
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,0);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);

% Define a fine grid of angle directions for 2-D search
varphi_range = linspace(-pi/3,pi/3,varphiSRes);
theta_range = linspace(-pi/3,pi/3,thetaSRes);

a_range = zeros(M,varphiSRes,thetaSRes); % [M,Azimuth,Elevation]
parfor i = 1:length(theta_range) % for each elevation
    a_range(:,:,i) = ...
        UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,varphiSRes),d_V,d_H);
end
load('widebeam16.mat');
widebeamresponses = [beamresponses,firsttarget,secondtarget];
%% Start the simulation
for n1 = 1:nbrOfAngleRealizations
    disp(n1);
    
    % generate channel vectors
    azimuth = unifrnd(-pi/3,pi/3,1);
    elevation = unifrnd(-pi/3,pi/3,1);
    g = UPA_Evaluate(lambda,M_V,M_H,azimuth,elevation,d_V,d_H);

    var_amp_d= 64;
    d = sqrt(var_amp_d/2) * (randn + 1i*randn);

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);
    for n2 = 1:nbrOfNoiseRealizations
        disp(n2)
        %Generate the noise
        noise = (randn(CBL,1)+1i*randn(CBL,1))/sqrt(2);

        %Select which two RIS configurations from the grid of beams to start with
        utilize = false(CBL,1);
        idx = unidrnd(CBL,2,1);
        while idx(1) == idx(2)
            idx(2) = unidrnd(CBL);
        end
        utilize(idx) = true;

        RISconfigs = Dh_angles*beamresponses(:,utilize);
        B = RISconfigs';

        %Go through iterations by adding extra RIS configurations in the estimation
        for itr = 1:Plim-1

            %Generate the received signal
            y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
 
            % Estimate the Channel using the developed MLE
            [~,var_phas_d_est,~,~,g_est,~,~] = MLE(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);

            %Estimate the RIS configuration that (approximately) maximizes the SNR
            RISconfig = angle(Dh*g_est)-var_phas_d_est;

            %Compute the corresponding achievable rate
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);

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
            end
        end
        SNR_total(:,n1,n2) = (abs(sqrt(SNR_pilot)*(B*Dh*g + d)).^2) ./ abs(noise(1:itr+1,1)).^2;
    
        %Select two wideband RIS configurations
        utilize = false(CBL+2,1);
        utilize(end-1:end) = true;

        RISconfigs = Dh_angles*widebeamresponses(:,utilize);
        B = RISconfigs';

        %Go through iterations by adding extra RIS configurations in the estimation
        for itr = 1:Plim-1

            %Generate the received signal
            y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
            
            % Estimate the Channel using the developed MLE
            [~,var_phas_d_est,~,~,g_est,~,~] = MLE(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);
    
            %Estimate the RIS configuration that (approximately) maximizes the SNR
            RISconfig = angle(Dh*g_est)-var_phas_d_est;

            %Compute the corresponding achievable rate
            rate_wide(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
            
            %Find an extra RIS configuration to use for pilot transmission
            if itr < Plim -1
                %Find which angles in the grid-of-beams haven't been used
                unusedBeamresponses = widebeamresponses;
                unusedBeamresponses(:,utilize==1) = 0;
                %Guess what the channel would be with the different beams
                guessBeam = Dh*unusedBeamresponses;

                %Find which of the guessed channels matches best with the
                %currently best RIS configuration
                closestBeam = abs(exp(-1i*RISconfig).'*guessBeam); 
                [~,bestUnusedBeamidx] = max(closestBeam);

                %Add a pilot transmission using the new RIS configuration
                utilize(bestUnusedBeamidx) = true;
                RISconfigs = Dh_angles*widebeamresponses(:,utilize);
                B = RISconfigs';
            end
        end
        SNR_w_total(:,n1,n2) = (abs(sqrt(SNR_pilot)*(B*Dh*g + d)).^2) ./ abs(noise(1:itr+1,1)).^2;
    end
end
rate = mean(mean(rate_proposed,3),2);
widerate = mean(mean(rate_wide,3),2);
plot(2:Plim,repelem(mean(capacity),1,Plim-1),'LineWidth',2);
hold on;
plot(2:Plim,widerate(1:Plim-1),'linewidth',2);
plot(2:Plim,rate(1:Plim-1),'LineWidth',2);
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
fig.Children.TickLabelInterpreter = 'latex';
legend('Capacity','Widebeam','NarrowBeam','Interpreter','latex');
grid on;

SNR_total_new = mean(mean(SNR_total,3),2);
SNR_w_total_new = mean(mean(SNR_w_total,3),2);
figure;
plot(pow2db(SNR_total_new),'LineWidth',2);
hold on;
plot(pow2db(SNR_w_total_new),'LineWidth',2);
legend('Narrowbeam','Widebeam');

save('draft2.mat','rate_proposed','rate_wide','SNR_total','SNR_w_total','capacity','Plim');