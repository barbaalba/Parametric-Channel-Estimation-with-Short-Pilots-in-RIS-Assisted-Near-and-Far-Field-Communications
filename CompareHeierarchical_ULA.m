clc;clear;close all;

%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
K = db2pow(12); % Rician K-factor in dB
Hierarchical_Search = true;
%ULA RIS configuration
M = 256;
d_H = 1/2; % in wavelnegth

if Hierarchical_Search
    W_Hierarch = ULA_Hierarchical(lambda,M,d_H);
end
%% Channel Estimation Parameters
% search resolution (It is very important)
varphiSRes = 4*M;
Plim = 2*log2(M); % number of pilots limited by Hierachical search

% Define a fine grid of angle directions and distance to be used in
% estimator
varphi_range = linspace(-pi/3,pi/3,varphiSRes);
% obtain the array response vectors for all azimuth
a_range = ULA_Evaluate(lambda,M,varphi_range,d_H);

%Set the SNR
SNRdB_pilot = 0;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -10;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;

% Generate channel between BS and RIS
h = ULA_Evaluate(lambda,M,varphi_BS,d_H);
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 500;
nbrOfNoiseRealizations = 10;

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_Hierchical = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
SNR = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
%% Iniitilize channel estimation 

% Generate code book
CBL = M;
beamresponses = fft(eye(M));
load('widebeamULA256.mat');
beamresponses = [W_Hierarch(:,:,log2(M)) W_Hierarch(:,1,1) W_Hierarch(:,2,1)];

for n1 = 1:nbrOfAngleRealizations
    disp(n1);

    azimuth = unifrnd(-pi/3,pi/3,1);

%     g = sqrt(K/(K+1)) * ULA_Evaluate(lambda,M,azimuth,d_H) + ...
%             sqrt(1/(K+1)/2)*(randn(M,1) + 1i*randn(M,1));
    g = ULA_Evaluate(lambda,M,azimuth,d_H);
            
    var_amp_d= 64;
    d = sqrt(var_amp_d/2) * (randn + 1i*randn);

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);

    for n2 = 1:nbrOfNoiseRealizations
        disp(n2);
        %Generate the noise
        noise = (randn(CBL,1)+1i*randn(CBL,1))/sqrt(2);
        utilize = false(CBL+2,1);
        utilize(end-1:end) = true;

        RISconfigs = Dh_angles*beamresponses(:,utilize);
        B = RISconfigs';   
        %Go through iterations by adding extra RIS configurations in the estimation
        for itr = 1:Plim-1

            %Generate the received signal
            y =  sqrt(SNR_pilot)*(B*Dh*g + d)  + noise(1:itr+1,1);
            SNR(itr,n1,n2) = mean(abs(y).^2);
            %[d_est,var_phas_d_est,var_amp_g_est,var_phas_g_est,g_est,Azidx,Elidx] = MLE(y,itr+1,B,Dh,a_range,SNR_pilot);
            [d_est,g_est] = MLE3DOptimized(y,itr+1,B,Dh,a_range,SNR_pilot);
                      
            RISconfig = exp(-1i * ( angle(Dh*g_est) - angle(d_est) ) );

            %Compute the corresponding achievable rate
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(RISconfig.'*Dh*g + d).^2);
     
            %Find an extra RIS configuration to use for pilot transmission
            if itr < Plim -1 

                %Find which angles in the grid-of-beams haven't been used
                unusedBeamresponses = beamresponses;
                unusedBeamresponses(:,utilize==1) = 0;

                %Guess what the channel would be with the different beams
                guessBeam = Dh_angles*unusedBeamresponses;

                %Find which of the guessed channels matches best with the
                %currently best RIS configuration
                closestBeam = abs(RISconfig.'*guessBeam); 
                [~,bestUnusedBeamidx] = max(closestBeam);

                %Add a pilot transmission using the new RIS configuration
                utilize(bestUnusedBeamidx) = true;
                RISconfigs = Dh_angles*beamresponses(:,utilize);
                B = RISconfigs';

            end
        end

        if Hierarchical_Search
            beamnum = 2;
            for k = 1:size(W_Hierarch,3)
                if k == 1
                    Beamidx = 1:beamnum;
                end
                B = (Dh_angles*W_Hierarch(:,Beamidx,k))';
                y = sqrt(SNR_pilot)*(B*Dh*g + d) + noise(beamnum*(k-1)+1:beamnum*k,1);
                [~,idx_max] = max(abs(y).^2);
                Bestidx = Beamidx(idx_max);
                RISconfig = B(idx_max,:);
                rate_Hierchical(k*beamnum,n1,n2) = log2(1+SNR_data*abs(RISconfig*Dh*g + d).^2);
%                 Beamidx(1:2) = 2*Bestidx-1:2*Bestidx; 
%                 Beamidx(3:4) = Beamidx(1:2) + 2^(k+1);
%                 Beamidx = Beamidx + 2^(k+1)*floor((Bestidx-1)/2^k);
                Beamidx = Beamidx(idx_max)*2-1:Beamidx(idx_max)*2;
            end
        end
        
    end
end

plot(2:Plim,repelem(mean(capacity),1,Plim-1),'--k','LineWidth',2);
hold on;
rate = mean(mean(rate_proposed,3),2);
rate_H = mean(mean(rate_Hierchical,3),2);
plot(2:Plim,rate(1:Plim-1),'LineWidth',2);
plot(beamnum:beamnum:Plim,rate_H(2:2:Plim));
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
fig.Children.TickLabelInterpreter = 'latex';
legend('Capacity','Widebeam beam initialization','Hierarchical Search','Interpreter','latex');
grid on; 
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);
