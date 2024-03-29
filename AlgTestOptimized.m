% This code is to test the functionality of the implemented algorithm to
% estimate the near field channel 
clc;clear;close all;
%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
NFConf = false; % True or false to specify which case to simulate (Near field or far field)
porpose = false; % True means it will estimate the channel through the proposed method
LSConf = false; % To include the LS estimation of the channel
Far_approx = true;
Racian = false; % The Rician channel model 
K = db2pow(8); % Rician K-factor in dB

%UPA Element configuration
M_H = 32; M_V = 32; M = M_H*M_V;
d_H = 1/2; d_V = 1/2; %In wavelengths
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D = sqrt(Hsize^2+Vsize^2); % Diagonal size of the array
d_fraun = 2*(D^2)/lambda; %2*diagonal^2/lambda
d_NF = d_fraun/10; % the upper threshold to be in near-field region
d_bjo = 2*D; % the lower threshold to be in near-field region with constant amplitude
disp(['The near-field upper threshold is ' num2str(d_NF) ' (m)']);
disp(['Bjornson distance is ' num2str(d_bjo) ' (m)']);
if d_NF < d_bjo
    disp('The code will not work since bjornson distance is greater than upper near field distance.');
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
R = UPAcorrelation(M_H,M_V,d_H,d_V,lambda); % correlation matrix

%% Channel Estimation Parameters
% search resolution (It is very important)
varphiSRes = 4*M_H;
thetaSRes = 2*M_V;
distSRes = M_H; % Distance resolution is very important to avoid convergence
Plim = 20; % number of pilots

%Set the SNR
SNRdB_pilot = -10;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -20;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H);
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 40;
nbrOfNoiseRealizations = 5;

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
Far_rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
SNR_generic = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
SNR_far = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
d_NMSE_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
g_NMSE_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
Far_d_NMSE_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
Far_g_NMSE_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
d_NMSE_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
v_NMSE_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
g_NMSE_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
%% Iniitilize channel estimation 
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,0);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);
% load and add two wide beam for 32 * 32 RIS with d_H = d_V = 1/2
load("WideTwobeam32.mat");
beamresponses = [beamresponses,firsttarget,secondtarget];
% use the same code book for far-field 
Farbeamresponses = beamresponses;
for n1 = 1:nbrOfAngleRealizations
    disp(n1);
    
    % generate the near-field channel 
    if NFConf
        d_t = unifrnd(d_bjo,d_NF,1);
    else
        d_t = unifrnd(d_NF,10*d_fraun);
    end
    azimuth = unifrnd(-pi/3,pi/3,1);
    elevation = unifrnd(-pi/3,pi/3,1);

    % Define a fine grid of angle directions and distance to be used in
    % estimator
    varphi_range = linspace(azimuth-pi/12,azimuth+pi/12,varphiSRes);
    theta_range = linspace(elevation-pi/24,elevation+pi/24,thetaSRes);
    dist_range = zeros(1,distSRes);
    if NFConf
        mind = max([d_bjo,d_t-d_bjo/8]);
        maxd = min([d_fraun,d_t+d_fraun/8]);
    else
        mind = max([d_NF,d_t-d_fraun/4]);
        maxd = min([10*d_fraun,d_t+d_fraun/4]);
    end
    dist_range(:) = linspace(mind,maxd,distSRes);

    if porpose
        % obtain the array response vectors for all azimuth-elevation-distance
        % triplet using the exact expression
        a_range = zeros(M,varphiSRes,thetaSRes,distSRes); % [M,Azimuth,Elevation,distance]
        for l =1:length(dist_range) % for each distance
            d_red = repelem(dist_range(l),1,varphiSRes);
            parfor i = 1:length(theta_range) % for each elevation
                a_range(:,:,i,l) = ...
                    nearFieldChan(d_red,varphi_range,repelem(theta_range(i),1,varphiSRes),U,lambda);
            end
        end
    end

    if Far_approx
        % Define a fine grid of angle directions to analyze when searching for
        % angle of arrival in far field
        a_FarAppx_range = zeros(M,varphiSRes,thetaSRes); % [M,Azimuth,Elevation]
        % obtain the array response vectors for all azimuth-elevation pairs
        parfor i = 1:thetaSRes
        a_FarAppx_range(:,:,i) = ...
        UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,varphiSRes),d_V,d_H);
        end
    end
    
    if ~Racian     
        g = nearFieldChan(d_t,azimuth,elevation,U,lambda); 
        var_amp_d= 64;
        d = sqrt(var_amp_d/2) * (randn + 1i*randn);
    else
        g = sqrt(K/(K+1)) * nearFieldChan(d_t,azimuth,elevation,U,lambda) + ...
            sqrtm(R)* sqrt(1/(K+1)/2)*(randn(M,1) + 1i*randn(M,1));
        var_amp_d= 64;
        d = sqrt(var_amp_d/2) * (randn + 1i*randn);
    end

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);

    for n2 = 1:nbrOfNoiseRealizations
        disp(n2);
        %Generate the noise
        noise = (randn(M+1,1)+1i*randn(M+1,1))/sqrt(2);
        if porpose
            % Estimate the channel using either all pilots or some of them
            %Select which two RIS configurations from the grid of beams to start with
            utilize = false(CBL+2,1);
            utilize(end-1:end) = true;
            
            RISconfigs = Dh_angles*beamresponses(:,utilize);
            B = RISconfigs';
    
            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1
    
                %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
                
                % Estimate the Channel using the developed MLE
                [d_est,var_phas_d_est,~,~,g_est,~,~] = MLE3D(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);
                
                % Estimation performance
                d_NMSE_proposed(itr,n1,n2) = norm(d_est - d)^2 / norm(d)^2;
                g_NMSE_proposed(itr,n1,n2) = norm(g_est - g)^2 / norm(g)^2;
    
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
            SNR_generic(:,n1,n2) = (abs(sqrt(SNR_pilot)*(B*Dh*g + d)).^2) ./ abs(noise(1:itr+1,1)).^2;
        end

        %% Far-Field approximation of the channel
        if Far_approx
            %Select which two RIS configurations from the grid of beams to start with
            utilize = false(CBL+2,1);
            utilize(end-1:end) = true;
            %Define the initial transmission setup
            RISconfigs = Dh_angles*Farbeamresponses(:,utilize);
            B = RISconfigs';
    
            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1
                %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
                
                % Estimate the Channel using the developed MLE
                [d_est,var_phas_d_est,~,~,g_est,~,~] = MLE(y,itr+1,B,Dh,a_FarAppx_range,SNR_pilot);
                
                % Estimation performance
                Far_d_NMSE_proposed(itr,n1,n2) = norm(d_est - d)^2 / norm(d)^2;
                Far_g_NMSE_proposed(itr,n1,n2) = norm(g_est - g)^2 / norm(g)^2;
             
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig = angle(Dh*g_est)-var_phas_d_est;
                
                %Compute the corresponding achievable rate
                Far_rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g + d)^2);
    
                %Find an extra RIS configuration to use for pilot transmission
                if itr < Plim-1
                    %Find which angles in the grid-of-beams haven't been used
                    unusedBeamresponses = Farbeamresponses;
                    unusedBeamresponses(:,utilize==1) = 0;
    
                    %Guess what the channel would be with the different beams
                    guessBeam = Dh*unusedBeamresponses;
    
                    %Find which of the guessed channels matches best with the
                    %currently best RIS configuration
                    closestBeam = abs(exp(-1i*RISconfig).'*guessBeam); 
                    [~,bestUnusedBeamidx] = max(closestBeam);
    
                    %Add a pilot transmission using the new RIS configuration
                    utilize(bestUnusedBeamidx) = true;
                    RISconfigs = Dh_angles*Farbeamresponses(:,utilize);
                    B = RISconfigs';
                end
            end
            SNR_far(:,n1,n2) = (abs(sqrt(SNR_pilot)*(B*Dh*g + d)).^2) ./ abs(noise(1:itr+1,1)).^2;
        end
        %% LS estimation
        if LSConf
            %LS estimation
            DFT = fft(eye(M+1));
            randomOrdering = randperm(M+1);
        
            %Go through iterations by adding extra RIS configurations in the estimation
            parfor itr = 1:Plim-1

                B_LS = transpose(DFT(:,randomOrdering(1:itr+1)));
    
                v = [d;Dh*g]; % cascaded channel and direct path 
                y = sqrt(SNR_pilot)*B_LS*v + noise(1:itr+1,1);
                
                %Compute LS estimate without parametrization
                v_LS = pinv(B_LS)*y/sqrt(SNR_pilot);
                % Estimation performance
                d_NMSE_LS(itr,n1,n2) = norm(v_LS(1) - d)^2 / norm(d)^2;
                v_NMSE_LS(itr,n1,n2) = norm(v_LS(2:end) - v(2:end))^2 / norm(v(2:end))^2;
                g_NMSE_LS(itr,n1,n2) = norm(inv(Dh)*v_LS(2:end) - g)^2/norm(g)^2
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig =  angle(v_LS(2:end)) - angle(v_LS(1));
    
                %Compute the corresponding achievable rate
                rate_LS(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig).'*Dh*g+d).^2);

            end
        end
    end
end

rate = mean(mean(rate_proposed,3),2);
rate_LS = mean(mean(rate_LS,3),2);
rate_FF = mean(mean(Far_rate_proposed,3),2);
plot(log2(2:Plim),repelem(mean(capacity),1,Plim-1),'LineWidth',2);
hold on;
plot(log2(2:Plim),rate(1:Plim-1),'LineWidth',2);
plot(log2(2:Plim),rate_FF(1:Plim-1),'LineWidth',2);
plot(log2(2:Plim),rate_LS(1:Plim-1),'LineWidth',2);
xlim([0 log2(Plim)]);
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
fig.Children.TickLabelInterpreter = 'latex';
legend('Capacity','Generic method','Far-Field approximation','Least Squares','Interpreter','latex');
grid on;
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
xticklabels({'$0$','$2^2$','$2^4$','$2^6$','$2^8$','$2^{10}$'});                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);


figure;
d_NMSE = mean(mean(d_NMSE_proposed,3),2);
g_NMSE = mean(mean(g_NMSE_proposed,3),2);
plot(2:Plim,pow2db(d_NMSE(1:Plim-1)),'LineWidth',2);
hold on;
plot(2:Plim,pow2db(g_NMSE(1:Plim-1)),'LineWidth',2);
plot(2:Plim,pow2db(mean(mean(d_NMSE_LS(1:Plim-1,:,:),3),2)));
plot(2:Plim,pow2db(mean(mean(g_NMSE_LS(1:Plim-1,:,:),3),2)));
ylabel('NMSE [dB]','FontSize',20,'Interpreter','latex');
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
legend('Direct channel $\mathbf{d}$','Channel to the RIS $\mathbf{g}$','Fontsize',20,'interpreter','latex');
grid on;
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);

save('FinalSISOLOS_LowSNR_NearFeild_together_new.mat','rate_proposed','rate_LS',...
      'Far_rate_proposed','capacity','Plim','SNR_generic','SNR_far',...
      'd_NMSE_proposed','g_NMSE_proposed','g_NMSE_LS','d_NMSE_LS');