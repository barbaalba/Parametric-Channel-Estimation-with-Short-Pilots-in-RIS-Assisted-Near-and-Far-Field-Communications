% 
clear; clc; close all;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

% RIS parameters
M_H = 16; M_V = 16; M = M_H*M_V;
M_H_red = 16; M_V_red = 4; M_red = M_H * M_V; % breaking section
d_H = 1/2;d_V = 1/2;

[ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V_red,M_H_red,d_V,d_H,0,0); 

% Consider below the RIS (Elevation <= 0)
elidx = find(ElAngles==0);
ElAngles = ElAngles(1:elidx);
fl = fieldnames(AzAngles);
for i = elidx+1:length(fl)
       AzAngles = rmfield(AzAngles,fl{i});
end

beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V_red,M_H_red,d_V,d_H);
samp = M_V/M_V_red;
sampl = M/samp;
firsttarget = zeros(M_V*M_H,1);
secondtarget = zeros(M_V*M_H,1);
for i = 1:samp 
    firsttarget((i-1)*sampl+1:i*sampl) = beamresponses(:,size(beamresponses,2)-i+1);
    secondtarget((i-1)*sampl+1:i*sampl) = beamresponses(:,i+1);
end

%% Plot 1
%Prepare to plot colors on a sphere
N = 540;
[X,Y,Z] = sphere(N);

%Prepare to compute channel gains on the sphere
gainMap = zeros(size(X));

%Go through all azimuth and elevation angles at specific distance
for n = 1:size(X,1)
    parfor m = 1:size(X,2)
        % To avoid plotting back of the antenna
        if X(n,m) < 0
            gainMap(n,m) = 0;
        else
            %Compute received power according to (7.28) in "Massive MIMO networks"
            [phi2,theta2] = cart2sph(X(n,m),Y(n,m),Z(n,m));
            gainMap(n,m) = abs(1/sqrt(M_H*M_V)*...
                firsttarget'*...
                UPA_Evaluate(lambda,M_V,M_H,phi2,theta2,d_V,d_H))^2;
        end
    end
end

% Plot beamforming gain
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
figure;
surf(X,Y,Z,pow2db(gainMap),'EdgeColor','none');
xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
caxis([0 pow2db(M_H*M_V)]); % limit the color map
colormap(flipud(hot));
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');
axis equal;
set(gca,'color',[0.9 0.9 0.9]);
set(gca,'fontsize',18);
view(122,30);
hold on;
xlim([0,1]); % To cut the sphere in half
varphiAngles = linspace(-pi/2,pi/2,100);
x_circ = cos(varphiAngles);
y_circ = sin(varphiAngles);
plot3(x_circ,y_circ,zeros(size(x_circ)),'k:','LineWidth',2);
%% plot 2 
%Prepare to plot colors on a sphere
N = 500;
[X,Y,Z] = sphere(N);

%Prepare to compute channel gains on the sphere
gainMap = zeros(size(X));

%Go through all azimuth and elevation angles at specific distance
for n = 1:size(X,1)
    parfor m = 1:size(X,2)
        % To avoid plotting back of the antenna
        if X(n,m) < 0
            gainMap(n,m) = 0;
        else
            %Compute received power according to (7.28) in "Massive MIMO networks"
            [phi2,theta2] = cart2sph(X(n,m),Y(n,m),Z(n,m));
            gainMap(n,m) = abs(1/sqrt(M_H*M_V)*...
                secondtarget'*...
                UPA_Evaluate(lambda,M_V,M_H,phi2,theta2,d_V,d_H))^2;  
        end
    end
end

% Plot beamforming gain
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
figure;
surf(X,Y,Z,pow2db(gainMap),'EdgeColor','none');
xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
caxis([0 pow2db(64)]); % limit the color map
colormap(flipud(hot));
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');
axis equal;
set(gca,'color',[0.9 0.9 0.9]);
set(gca,'fontsize',18);
view(122,30);
hold on;
xlim([0,1]); % To cut the sphere in half

%%

firsttarget = beamresponses(:,end-samp+1:end);
firsttarget = sum(firsttarget,2);
firsttarget = firsttarget / norm(firsttarget) * sqrt(M_V*M_H);
%Prepare to compute channel gains on the sphere
gainMap = zeros(size(X));

%Go through all azimuth and elevation angles at specific distance
for n = 1:size(X,1)
    parfor m = 1:size(X,2)
        % To avoid plotting back of the antenna
        if X(n,m) < 0
            gainMap(n,m) = 0;
        else
            %Compute received power according to (7.28) in "Massive MIMO networks"
            [phi2,theta2] = cart2sph(X(n,m),Y(n,m),Z(n,m));
            gainMap(n,m) = abs(1/sqrt(M_H*M_V)*...
                firsttarget'*...
                UPA_Evaluate(lambda,M_V,M_H,phi2,theta2,d_V,d_H))^2;  
        end
    end
end

% Plot beamforming gain
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
figure;
surf(X,Y,Z,pow2db(gainMap),'EdgeColor','none');
xlabel('$x$','Interpreter','Latex');
ylabel('$y$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
caxis([0 pow2db(M_H*M_V)]); % limit the color map
colormap(flipud(hot));
hBar = colorbar;
set(hBar, 'TickLabelInterpreter', 'latex');
axis equal;
set(gca,'color',[0.9 0.9 0.9]);
set(gca,'fontsize',18);
view(122,30);
hold on;
xlim([0,1]); % To cut the sphere in half
% Plot the dot lines around the sphere
