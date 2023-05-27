close all;clear;clc;rng(9);
%% Room Scenario Initialization
Xmax = 2; % Room size = [-Xmax, Xmax]
Ymax = 2; % Room size = [-Ymax, Ymax]
numUE = 1; % number of users
RWL = 2;  % Length of the random walk
Speedlim = [0.1,0.4]; % meter per second movement
RIS_center = [0,0,0]; % The RIS center point

%% Rx/Tx/RIS Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

% RIS parameters
M_H = 32; M_V = 32; M = M_H * M_V;
d_H = 1/2;d_V = 1/2;
% RIS element position
i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% The normal convention to evaluate position of the elements in the room
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [RIS_center(1); RIS_center(2)+ym; RIS_center(3)+zm]; % Position of the m-th element
end
% Near/Far-field parameters
Hsize = M_H*d_H*lambda;
Vsize = M_V*d_V*lambda;
D = sqrt(Hsize^2+Vsize^2); % Diagonal size of the array
d_fraun = 2*(D^2)/lambda; %2*diagonal^2/lambda
d_NF = d_fraun/10; % the upper threshold to be in near-field region
d_bjo = 2*D; % the lower threshold to be in near-field region with constant amplitude
disp(['Lower and upper distance are: ' num2str(d_bjo) ' , '  num2str(d_NF) ' (m)']);

%% Random walk 
plt = true; % To plot the trajectory
pltconf = 'continous'; % 'continous' or 'discrete'
[x_t,y_t] = randomwalk(numUE,RWL,Xmax,Ymax,Speedlim,d_bjo,RIS_center);
z_t = repelem(-0.1,numUE,size(x_t,2)); % Z-coordinate of the user(does not change)
[azimuth,elevation,Cph,d_t] = ChanParGen(x_t,y_t,z_t,RIS_center,lambda);
if plt
    plotTrajectory(x_t,y_t,azimuth,elevation,pltconf,Xmax,Ymax,RIS_center,d_bjo,d_NF);
end
% near-field percentage
NF_percent = sum(d_t < d_NF & d_t > d_bjo,'all')/size(d_t,2)/numUE*100;
disp(['How much percentage in near field: ' num2str(NF_percent)]);

%% Real Channel generator

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Generate channels
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H); % [M,1]
Dh = diag(h);
Dh_angles = diag(h./abs(h));

g_t = realChanGen(x_t,y_t,z_t,U,lambda); % [M,T]

%% Initialize the generic channel estimator 
nbrOfNoiseRealizations = 1;

i = @(m) mod(m-1,M_H); % Horizontal index
j = @(m) floor((m-1)/M_H); % Vertical index
U = zeros(3,M); % Matrix containing the position of the RIS antennas 
% The normal convention to evaluate position of the elements 
for m = 1:M  
    ym = (-(M_H-1)/2 + i(m))*d_H*lambda; % dislocation with respect to center in y direction
    zm = (-(M_V-1)/2 + j(m))*d_V*lambda; % dislocation with respect to center in z direction
    U(:,m) = [0; ym; zm]; % Relative position of the m-th element with respect to the center
end
%Set the SNR
SNRdB_pilot = 10;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = 0;
SNR_data = db2pow(SNRdB_data);
Plim = M_H/2; % number of pilots

% Channel estimation codebook
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,1.3264);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);

% grid search resolution (It is very important)
varphiSRes = 2*M_H;
thetaSRes = 2*M_V;
distSRes = 4*M_H; % Distance resolution is very important to avoid convergence

% Define a fine grid of angle directions and distance to be used in
% estimator
varphi_range = linspace(-pi/2,pi/2,varphiSRes);
theta_range = linspace(-pi/2,pi/2,thetaSRes);
dist_range = zeros(1,distSRes);
dist_range(1:distSRes/2) = linspace(d_bjo,d_NF,distSRes/2);
dist_range(distSRes/2+1:end) = linspace(d_fraun/9,10*d_fraun,distSRes/2);

a_range = zeros(M,varphiSRes,thetaSRes,distSRes); % [M,Azimuth,Elevation,distance]
% obtain the array response vectors for all azimuth-elevation-distance
% triplet using the exact expression
for l =1:length(dist_range) % for each distance
    d_t = repelem(dist_range(l),1,varphiSRes);
    parfor i = 1:length(theta_range) % for each elevation
        a_range(:,:,i,l) = ...
            nearFieldChan(d_t,varphi_range,repelem(theta_range(i),1,varphiSRes),U,lambda);
    end
end
%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,size(x_t,2));
rate_proposed = zeros(Plim,size(x_t,2),nbrOfNoiseRealizations);

%% Perform channel estimation
for n1=1:size(x_t,2)
    disp(n1);
    % current channel to be estimated
    g = g_t(:,n1);
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

%% Plot the result
rate = mean(mean(rate_proposed,3),2);
plot(2:Plim,repelem(mean(capacity),1,Plim-1),'LineWidth',2);
hold on;
plot(2:Plim,rate(1:Plim-1),'LineWidth',2);
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
fig.Children.TickLabelInterpreter = 'latex';
legend('Capacity','MLE','Interpreter','latex');
grid on;