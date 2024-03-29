clc;clear;close all;

%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
Hierarchical_Search = true;
%UPA Element configuration
M_H = 16; M_V = 16; M = M_H*M_V;
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

if Hierarchical_Search
    W_Hierarch = UPA_Hierarchical(M_H,M_V,d_H,d_V);
end
%% Channel Estimation Parameters
% search resolution (It is very important)
varphiSRes = 16*M_H;
thetaSRes = 2*M_H; % We consider far-field and elevation 0 
Plim = 4*log2(M_H); % number of pilots limited by Hierachical search

% Define a fine grid of angle directions and distance to be used in
% estimator
varphi_range = linspace(-pi/3,pi/3,varphiSRes);
theta_range = linspace(-pi/3,0,thetaSRes);
% obtain the array response vectors for all azimuth-elevation-distance
% triplet using the exact expression
a_range = zeros(M,varphiSRes,thetaSRes); % [M,Azimuth,Elevation]
parfor i = 1:length(theta_range) % for each elevation
    a_range(:,:,i) = ...
        UPA_Evaluate(lambda,M_V,M_H,varphi_range,repelem(theta_range(i),1,varphiSRes),d_V,d_H);
end

%Set the SNR
SNRdB_pilot = 0;
SNR_pilot = db2pow(SNRdB_pilot);

SNRdB_data = -10;
SNR_data = db2pow(SNRdB_data);

%Select angle to the base station (known value)
varphi_BS = -pi/6;
theta_BS = 0;

% Generate channel between BS and RIS
h = UPA_Evaluate(lambda,M_V,M_H,varphi_BS,theta_BS,d_V,d_H); % for signle BS
Dh = diag(h);
Dh_angles = diag(h./abs(h));

nbrOfAngleRealizations = 100;
nbrOfNoiseRealizations = 5;

%Save the rates achieved at different iterations of the algorithm
capacity = zeros(1,nbrOfAngleRealizations);
capacity_Hierchical = zeros(1,nbrOfAngleRealizations);
rate_proposed = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
rate_Hierchical = zeros(Plim,nbrOfAngleRealizations,nbrOfNoiseRealizations);
%% Iniitilize channel estimation 

% Generate code book
[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,pi/2,0);
beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);

load("WideTwobeam16.mat");
widebeamresponses = conj(Dh_angles * [beamresponses,firsttarget,secondtarget]);

for n1 = 1:nbrOfAngleRealizations
    disp(n1);

    azimuth = unifrnd(-pi/3,pi/3,1);
    elevation = unifrnd(-pi/3,0,1);

    g = UPA_Evaluate(lambda,M_V,M_H,azimuth,elevation,d_V,d_H);
            
    %Compute the exact capacity for the system Eq. (3)
    capacity(n1) = log2(1+SNR_data*sum(abs(Dh*g))^2);

    for n2 = 1:nbrOfNoiseRealizations
        disp(n2);
        %Generate the noise
        noise = (randn(CBL,1)+1i*randn(CBL,1))/sqrt(2);

        utilize = false(CBL+2,1);
        utilize(end-1:end) = true;

        RISconfigs = widebeamresponses(:,utilize);
        B = RISconfigs.';   
        %Go through iterations by adding extra RIS configurations in the estimation
        for itr = 1:Plim-1

            %Generate the received signal
            y =  sqrt(SNR_pilot)*B*Dh*g  + noise(1:itr+1,1);
            
            utilityfunction = zeros(varphiSRes,thetaSRes); % The utility function
            parfor i = 1:thetaSRes
                utilityfunction(:,i) = abs(y'*B*Dh*a_range(:,:,i)).^2./sum(abs(B*Dh*a_range(:,:,i)).^2,1);
            end
            %Extract the angle estimate
            [~,maxind] = max(utilityfunction,[],'all');
            [Azidx,Elidx] = ind2sub([varphiSRes,thetaSRes],maxind);
            a = a_range(:,Azidx,Elidx);
            
            RISconfig = exp( -1i * angle(Dh*a)  );

            %Compute the corresponding achievable rate
            rate_proposed(itr,n1,n2) = log2(1+SNR_data*abs(RISconfig.'*Dh*g).^2);
     
            %Find an extra RIS configuration to use for pilot transmission
            if itr < Plim -1 

                %Find which angles in the grid-of-beams haven't been used
                unusedBeamresponses = widebeamresponses;
                unusedBeamresponses(:,utilize==1) = 0;
                %Guess what the channel would be with the different beams
                guessBeam = unusedBeamresponses;

                %Find which of the guessed channels matches best with the
                %currently best RIS configuration
                closestBeam = abs(RISconfig'*guessBeam); 
                [~,bestUnusedBeamidx] = max(closestBeam);

                %Add a pilot transmission using the new RIS configuration
                utilize(bestUnusedBeamidx) = true;
                RISconfigs = widebeamresponses(:,utilize);
                B = RISconfigs.';

            end
        end

        if Hierarchical_Search
            for k = 1:size(W_Hierarch,3)
                if k == 1
                    Beamidx = 1:4;
                end
                B = (Dh_angles*W_Hierarch(:,Beamidx,k))';
                y = sqrt(SNR_pilot)*B*Dh*g + noise(4*(k-1)+1:4*k,1);
                [~,idx_max] = max(abs(y).^2);
                Bestidx = Beamidx(idx_max);
                RISconfig = B(idx_max,:);
                rate_Hierchical(k*4,n1,n2) = log2(1+SNR_data*abs(RISconfig*Dh*g).^2);
                Beamidx(1:2) = 2*Bestidx-1:2*Bestidx; 
                Beamidx(3:4) = Beamidx(1:2) + 2^(k+1);
                Beamidx = Beamidx + 2^(k+1)*floor((Bestidx-1)/2^k);
                %Beamidx = Beamidx(idx_max)*2-1:Beamidx(idx_max)*2;
            end
        end
        
    end
end
save('HierarchicalUPA_withoutDirect.mat','rate_proposed','capacity','rate_Hierchical')
plot(2:Plim,repelem(mean(capacity),1,Plim-1),'--k','LineWidth',2);
hold on;
rate = mean(mean(rate_proposed,3),2);
rate_H = mean(mean(rate_Hierchical,3),2);
plot(2:Plim,rate(1:Plim-1),'LineWidth',2);
plot(4:4:Plim,rate_H(4:4:Plim));
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
