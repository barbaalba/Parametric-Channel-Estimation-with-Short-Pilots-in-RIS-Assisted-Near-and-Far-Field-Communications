clc;clear;close all;

%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
K = db2pow(8); % Rician K-factor in dB
LS = true; % Least Squares Estiamtor
proposed = false; % The proposed MLE method

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
% BS config
N = 64; d_H_BS = 1/2;

%% Channel Estimation Parameters
% search resolution (It is very important)
varphiSRes = 16*M_H;
thetaSRes = 2*M_V;
distSRes = 1; % Distance resolution is very important to avoid convergence
Plim = 100; % number of pilots

%Set the SNR
SNRdB_pilot = -20;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -30;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Select the angle from the BS (Known value)
varphi_RIS = pi/4;

% Generate channel between BS and RIS
H = sqrt(K/(K+1)) * ULA_Evaluate(N,varphi_RIS,d_H_BS) * UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H)' +...
    sqrt(1/(K+1)/2)*(randn(N,M) + 1i*randn(N,M));
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H); % for signle BS
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 5;
nbrOfNoiseRealizations = 5;

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
d_NMSE_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
g_NMSE_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
d_NMSE_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
g_NMSE_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_LS = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
%% Iniitilize channel estimation 

% Generate code book
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,0);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);
if proposed
    load("WideTwobeam32.mat");
    widebeamresponses = conj([beamresponses,firsttarget,secondtarget]);
end

for n1 = 1:nbrOfAngleRealizations
    disp(n1);

    d_t = unifrnd(d_fraun,10*d_fraun);
    azimuth = unifrnd(-pi/3,pi/3,1);
    elevation = unifrnd(-pi/3,pi/3,1);
    if proposed
        % Define a fine grid of angle directions and distance to be used in
        % estimator
        varphi_range = linspace(azimuth-pi/6,azimuth+pi/6,varphiSRes);
        theta_range = linspace(elevation-pi/24,elevation+pi/24,thetaSRes);
        dist_range = zeros(1,distSRes);
        mind = max([d_bjo,d_t-d_bjo/8]);
        maxd = min([10*d_fraun,d_t+d_fraun/8]);
        dist_range(:) = d_t;%linspace(mind,maxd,distSRes);
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
    g = sqrt(K/(K+1)) * nearFieldChan(d_t,azimuth,elevation,U,lambda) + ...
            sqrtm(R)* sqrt(1/(K+1)/2)*(randn(M,1) + 1i*randn(M,1));
    var_amp_d= 64;
    d = sqrt(var_amp_d) * (sqrt(K/(K+1)) * ULA_Evaluate(N,unifrnd(-pi/3,pi/3,1),d_H_BS) + ...
        sqrt(1/2/(K+1))*(randn(N,1) + 1i*randn(N,1)));
    %d = sqrt(var_amp_d) * (randn(N,1) + 1i*randn(N,1));

    optRIS = altopt(H,g,d);
    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*norm(H*diag(g)*optRIS + d)^2);

    for n2 = 1:nbrOfNoiseRealizations
        disp(n2);
        %Generate the noise
        noise = (randn(CBL*N,1)+1i*randn(CBL*N,1))/sqrt(2);
        if proposed
            utilize = false(CBL+2,1);
            utilize(end-1:end) = true;
%             idx = randi(CBL,2,1);
%             while idx(2) == idx(1)
%                 idx(2) = randi(CBL,1);
%             end
%             utilize(idx) = true;
    
            RISconfigs = Dh_angles*widebeamresponses(:,utilize);
            B = RISconfigs.';   
            %Go through iterations by adding extra RIS configurations in the estimation
            for itr = 1:Plim-1

                F = kron(ones(itr+1,1),H) .* kron(B,ones(N,1)); 
                %Generate the received signal
                y =  sqrt(SNR_pilot) * ( F*g + kron(ones(itr+1,1),d) ) + noise(1:(itr+1)*N,1);
                
                % Estimate the Channel using the developed MLE
                [d_est,g_est]  = MLE4D(y,itr+1,F,a_range,SNR_pilot);
            
                % Estimation performance
                d_NMSE_proposed(itr,n1,n2) = norm(d_est - d)^2 / norm(d)^2;
                g_NMSE_proposed(itr,n1,n2) = norm(g_est - g)^2 / norm(g)^2;
    
                %Estimate the RIS configuration that (approximately) maximizes the SNR
                RISconfig = altopt(H,g_est,d_est);
    
                %Compute the corresponding achievable rate
                rate_proposed(itr,n1,n2) = log2(1+SNR_data*norm(H*diag(g)*RISconfig + d)^2);
         
                %Find an extra RIS configuration to use for pilot transmission
                if itr < Plim -1 

                    %Find which angles in the grid-of-beams haven't been used
                    unusedBeamresponses = widebeamresponses;
                    unusedBeamresponses(:,utilize==1) = 0;
                    %Guess what the channel would be with the different beams
                    guessBeam = Dh*unusedBeamresponses;

                    %Find which of the guessed channels matches best with the
                    %currently best RIS configuration
                    closestBeam = abs(RISconfig'*guessBeam); 
                    [~,bestUnusedBeamidx] = max(closestBeam);


                    %Add a pilot transmission using the new RIS configuration
                    utilize(bestUnusedBeamidx) = true;
                    RISconfigs = Dh_angles*widebeamresponses(:,utilize);
                    B = RISconfigs.';

                end
            end
        end

        if LS
            %LS estimation
%             DFT = DFTBookBuild(M_H,M_V);
            randomOrdering = randperm(CBL);
            chan = zeros(M+N,1);
            chan(1:N) = d;
            chan(N+1:end) = g;
           
            parfor itr = 1:Plim-1
                %B = transpose(DFT(:,randomOrdering(1:itr+1)));
                B = transpose(Dh_angles*conj(beamresponses(:,randomOrdering(1:itr+1))));
                F = kron(ones(itr+1,1),H) .* kron(B,ones(N,1));
                
                y =  sqrt(SNR_pilot) * ( F*g + kron(ones(itr+1,1),d) ) + noise(1:(itr+1)*N,1);
                Fnew = zeros(N*(itr+1),M+N);
                Fnew(:,1:N) = kron(ones(itr+1,1),eye(N));
                Fnew(:,N+1:end) = F;
               % y = sqrt(SNR_pilot)*Fnew*chan + noise(1:(itr+1)*N,1);
                chanest = pinv(Fnew)*y/sqrt(SNR_pilot);
    
                d_NMSE_LS(itr,n1,n2) = norm(chanest(1:N) - d)^2/norm(d)^2;
                g_NMSE_LS(itr,n1,n2) = norm(chanest(N+1:end) - g)^2/norm(g)^2;
    
                RISconfig_LS = altopt(H,chanest(N+1:end),chanest(1:N));
                %Compute the corresponding achievable rate
                rate_LS(itr,n1,n2) = log2(1+SNR_data*norm(H*diag(g)*RISconfig_LS + d)^2);
            end
        end
    end
end
%save('SingleAntennaFarField_Rician8_LowSNR_2.mat','Plim','rate_proposed','capacity','d_NMSE_proposed','g_NMSE_proposed');
plot(2:Plim,repelem(mean(capacity),1,Plim-1),'--k','LineWidth',2);
rate = mean(mean(rate_proposed,3),2);
hold on;
plot(2:Plim,rate(1:Plim-1),'LineWidth',2);
rate = mean(mean(rate_LS,3),2);
plot(2:Plim,rate(1:Plim-1),'LineWidth',2);
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
ylabel('Spectral Efficiency [b/s/Hz/]','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
fig.Children.TickLabelInterpreter = 'latex';
legend('Capacity','Widebeam beam initialization','LS','Interpreter','latex');
grid on; 
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);

figure;
d_NMSE = mean(mean(d_NMSE_proposed,3),2);
g_NMSE = mean(mean(g_NMSE_proposed,3),2);
plot(2:Plim,pow2db(d_NMSE(1:Plim-1)),'LineWidth',2);
hold on;
plot(2:Plim,pow2db(g_NMSE(1:Plim-1)),'LineWidth',2);
d_NMSE = mean(mean(d_NMSE_LS,3),2);
g_NMSE = mean(mean(g_NMSE_LS,3),2);
plot(2:Plim,pow2db(d_NMSE(1:Plim-1)),'LineWidth',2);
plot(2:Plim,pow2db(g_NMSE(1:Plim-1)),'LineWidth',2);
ylabel('NMSE [dB]','FontSize',20,'Interpreter','latex');
xlabel('Number of pilots','FontSize',20,'Interpreter','latex');
fig = gcf;
fig.Children.FontSize = 20;
legend('Direct channel $\mathbf{d}$','Channel to the RIS $\mathbf{g}$','Direct channel $\mathbf{d}$ - LS','Channel to the RIS $\mathbf{g}$ - LS','Fontsize',20,'interpreter','latex');
grid on;
ax = gca; % to get the axis handle
ax.XLabel.Units = 'normalized'; % Normalized unit instead of 'Data' unit 
ax.Position = [0.15 0.15 0.8 0.8]; % Set the position of inner axis with respect to
                           % the figure border
ax.XLabel.Position = [0.5 -0.07]; % position of the label with respect to 
                                  % axis
fig = gcf;
set(fig,'position',[60 50 900 600]);