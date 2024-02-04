close all;clear;clc;rng(10);
%% Room Scenario Initialization
Xmax = 4; % Room size = [0, Xmax]
Ymax = 2; % Room size = [-Ymax, Ymax]
numUE = 1; % number of users
RWL = 3;  % Length of the random walk in second
brkfreq = 1000; % number of channel instance per second
Speedlim = [0.9,1.1]; % meter per second movement
RIS_center = [0,0,0]; % The RIS center point
fest = brkfreq/10; % Estimation frequency
ftrack = brkfreq/100; % tracking frequency 
%% Rx/Tx/RIS Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

% RIS parameters
M_H = 32; M_V = 32; M = M_H * M_V;
d_H = 1/2;d_V = 1/2;
% Near/Far-field parameters
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D = sqrt(Hsize^2+Vsize^2); % Diagonal size of the array
d_fraun = 2*(D^2)/lambda; %2*diagonal^2/lambda
d_NF = d_fraun/10; % the upper threshold to be in near-field region
d_bjo = 2*D; % the lower threshold to be in near-field region with constant amplitude
disp(['Lower and upper distance are: ' num2str(d_bjo) ' , '  num2str(d_NF) ' (m)']);
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% The normal convention to evaluate position of the elements 
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end
%% Random walk 
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Speedlim,d_bjo,d_NF,RIS_center,brkfreq);
z_t = repelem(0,numUE,size(x_t,2)); % Z-coordinate of the user(does not change)
[azimuth,elevation,Cph,d_t] = ChanParGen(x_t,y_t,z_t,RIS_center,lambda);

% near-field percentage
NF_percent = sum(d_t < d_NF & d_t > d_bjo,'all')/size(d_t,2)/numUE*100;
disp(['How much percentage in near field: ' num2str(NF_percent)]);

%% Channel estimation parameter
% grid search resolution (It is very important)
varphiSRes = M_H;
thetaSRes = M_H/2;
distSRes = M_H; % Distance resolution is very important to avoid convergence
Plim = 10; % number of pilots
ptrack = 6; % number of tracking pilots

%Set the SNR
SNRdB_pilot = -10;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -20;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H); % [M,1]
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = size(x_t,2);
nbrOfNoiseRealizations = 5;

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
RISconfig = zeros(M,nbrOfNoiseRealizations);
%% Iniitilize channel estimation  
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,0);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);
trackresponses = beamresponses;
% load and add two wide beam for 32 * 32 RIS with d_H = d_V = 1/2
load("WideTwobeam32.mat");
beamresponses = [beamresponses,firsttarget,secondtarget];

%% Perform channel estimation
for n1=1:nbrOfAngleRealizations
    disp(n1);
    
    % current channel to be estimated
    g = nearFieldChan(d_t(n1),azimuth(n1),elevation(n1),U,lambda);
    var_amp_d= 64;
    d = sqrt(var_amp_d/2) * (randn + 1i*randn);

    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*(sum(abs(Dh*g)) + abs(d)).^2);

    for n2 = 1:nbrOfNoiseRealizations
        disp(n2);
        % Estimate the channel using either all pilots or some of them
        %Select which two RIS configurations from the grid of beams to start with
        if  n1 == 1

            % Define a fine grid of angle directions and distance to be used in
            % estimator
            varphi_range = linspace(max(azimuth(n1)-pi/16,-pi/2),min(azimuth(n1)+pi/16,pi/2),varphiSRes);
            theta_range = linspace(-pi/18,pi/18,thetaSRes);
            dist_range = zeros(1,distSRes);
            dist_range(:) = linspace(d_bjo,sqrt(Xmax^2+Ymax^2),distSRes);
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

            utilize = false(CBL+2,1);
            utilize(end-1:end) = true;

            RISconfigs = Dh_angles*beamresponses(:,utilize);
            B = RISconfigs'; %[P,M]

            %Generate the noise
            noise = (randn(CBL,1)+1i*randn(CBL,1))/sqrt(2);
    
            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1
    
                %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
                
                % Estimate the Channel using the developed MLE
                [~,var_phas_d_est,~,~,g_est,~,~] = MLE3D(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);
                
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig(:,n2) = angle(Dh*g_est)-var_phas_d_est;
    
                %Compute the corresponding achievable rate
                rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig(:,n2)).'*Dh*g + d)^2);
         
                %Find an extra RIS configuration to use for pilot transmission
                if itr < Plim 

                    %Find which angles in the grid-of-beams haven't been used
                    unusedBeamresponses = beamresponses;
                    unusedBeamresponses(:,utilize==1) = 0;
                    %Guess what the channel would be with the different beams
                    guessBeam = Dh*unusedBeamresponses;
    
                    %Find which of the guessed channels matches best with the
                    %currently best RIS configuration
                    closestBeam = abs(exp(-1i*RISconfig(:,n2)).'*guessBeam); 
                    [~,bestUnusedBeamidx] = max(closestBeam);
    
                    %Add a pilot transmission using the new RIS configuration
                    utilize(bestUnusedBeamidx) = true;
                    RISconfigs = Dh_angles*beamresponses(:,utilize);
                    B = RISconfigs';  

                end
            end
        elseif mod(n1,ftrack) == 1

            % Define a fine grid of angle directions and distance to be used in
            % estimator
            varphi_range = linspace(max(azimuth(n1)-pi/16,-pi/2),min(azimuth(n1)+pi/16,pi/2),varphiSRes);
            theta_range = linspace(-pi/18,pi/18,thetaSRes);
            dist_range = zeros(1,distSRes);
            dist_range(:) = linspace(d_bjo,sqrt(Xmax^2+Ymax^2),distSRes);
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
            val = abs(exp(-1i*RISconfig(:,n2)).'*Dh_angles*trackresponses);
            [~,idx] = sort(val,'descend');
            utilize = false(CBL,1);
            utilize(idx(1:2)) = true;
    
            RISconfigs = Dh_angles*trackresponses(:,utilize);
            B = RISconfigs'; %[P,M]

            %Generate the noise
            noise = (randn(CBL,1)+1i*randn(CBL,1))/sqrt(2);

            for itr = 1:ptrack-1
                %Generate the received signal
                y =  sqrt(SNR_pilot)*(B*Dh*g + d) + noise(1:itr+1,1);
                    
                % Estimate the Channel using the developed MLE
                [~,var_phas_d_est,~,~,g_est,~,~] = MLE3D(y,itr+1,B,Dh,a_range,lambda,M_V,M_H,d_V,d_H,varphi_range,theta_range,SNR_pilot);
                    
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig(:,n2) = angle(Dh*g_est)-var_phas_d_est;
        
                %Compute the corresponding achievable rate
                rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig(:,n2)).'*Dh*g + d)^2);
     
                %Find an extra RIS configuration to use for pilot transmission
                if itr < ptrack-1 
        
                    %Find which angles in the grid-of-beams haven't been used
                    unusedBeamresponses = trackresponses;
                    unusedBeamresponses(:,utilize==1) = 0;
                    %Guess what the channel would be with the different beams
                    guessBeam = Dh*unusedBeamresponses;
        
                    %Find which of the guessed channels matches best with the
                    %currently best RIS configuration
                    closestBeam = abs(exp(-1i*RISconfig(:,n2)).'*guessBeam); 
                    [~,bestUnusedBeamidx] = max(closestBeam);
        
        
                    %Add a pilot transmission using the new RIS configuration
                    utilize(bestUnusedBeamidx) = true;
                    RISconfigs = Dh_angles*trackresponses(:,utilize);
                    B = RISconfigs';
        
                else
                    rate_proposed(Plim-1,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig(:,n2)).'*Dh*g + d)^2);
                end
            end
        else
             %Compute the corresponding achievable rate
             rate_proposed(Plim-1,n1,n2) = log2(1+SNR_data*abs(exp(-1i*RISconfig(:,n2)).'*Dh*g + d)^2);
        end

    end
end
save('FinalRandomWalk0est10trck4p.mat','rate_proposed','capacity','Plim','ptrack',...
    'fest','ftrack','brkfreq','RWL','azimuth','d_t','nbrOfAngleRealizations');
%% Plot the result
rate = mean(rate_proposed(Plim-1,:,:),3);
figure;
plot(1:nbrOfAngleRealizations,capacity,'LineWidth',2,'LineStyle','--','Color','k');
hold on;
plot(1:nbrOfAngleRealizations,rate,'LineWidth',2);
xticks(0:brkfreq:RWL*brkfreq);
yticks(0:2:20);
xlim([0 RWL*brkfreq]);
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'});
% plot(1:fest:nbrOfAngleRealizations,rate(1:fest:nbrOfAngleRealizations),...
%     'LineStyle','none','Marker','o','MarkerSize',10,'LineWidth',2);
ax = gca;
ax.TickLabelInterpreter = 'latex';
fig = gcf;
fig.Children.FontSize = 20;
xlabel('Time [s]','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
legend('Capacity','Live Performace','Re-estimation','FontSize',20,'Interpreter','latex');
grid on;
ylim([0,18]);
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]); %[left bottom width height]
