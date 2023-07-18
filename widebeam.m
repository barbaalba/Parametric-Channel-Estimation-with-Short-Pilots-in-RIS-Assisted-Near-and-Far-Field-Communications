clear; clc; close all;
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength

% RIS parameters
M_H = 16; M_V = 16; M = M_H*M_V;
d_H = 1/2;d_V = 1/2;
conf  = 1; % 1= widebam, 2= Hedgehog shape beam
if conf == 1   
    N_s = sqrt(2*M/(pi*d_V*d_H));
    N_s = 2^floor(log2(N_s)); % number of elements per sector
    ax  = 1*(M_H<=M_V);
    % if horizental axis is smaller
    if ax 
        N_h = M_H; N_v = N_s/N_h; % structure of each sector
    % if vertical axis is smaller
    else
        N_v = M_V;N_h = N_s/N_v;
    end
    s = M/N_s; % number of sectors
    [ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_H,M_V,d_V,d_H,0,0); 
    
    beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);
    firsttarget = zeros(M,1);
    secondtarget = zeros(M,1);
    % per section config association
    for i = 1:s 
        firsttarget((i-1)*N_s+1:i*N_s) = beamresponses((i-1)*N_s+1:i*N_s,size(beamresponses,2)-i+1);
        secondtarget((i-1)*N_s+1:i*N_s) = beamresponses((i-1)*N_s+1:i*N_s,i);
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
    caxis([0 pow2db(N_s)]); % limit the color map
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
    caxis([0 pow2db(N_s)]); % limit the color map
    colormap(flipud(hot));
    hBar = colorbar;
    set(hBar, 'TickLabelInterpreter', 'latex');
    axis equal;
    set(gca,'color',[0.9 0.9 0.9]);
    set(gca,'fontsize',18);
    view(122,30);
    hold on;
    xlim([0,1]); % To cut the sphere in half
    x_circ = cos(varphiAngles);
    y_circ = sin(varphiAngles);
    plot3(x_circ,y_circ,zeros(size(x_circ)),'k:','LineWidth',2);
    save('widebeam16.mat','firsttarget','secondtarget');
else 
    [ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,0,0); 
    % Consider below the RIS (Elevation <= 0)
    elidx = find(ElAngles==0);
    ElAngles = ElAngles(1:elidx);
    fl = fieldnames(AzAngles);
    for i = elidx+1:length(fl)
           AzAngles = rmfield(AzAngles,fl{i});
    end
    fl = fieldnames(AzAngles);
    Azright = struct;
    for i = 1:length(fl)
        azidx = find(AzAngles.(fl{i})==0);
        Azright.(fl{i}) = (AzAngles.(fl{i})(azidx:end));
    end
    rbeamresponses = UPA_Codebook(lambda,ElAngles,Azright,M_V,M_H,d_V,d_H);
    rsect = 2^floor(log2(size(rbeamresponses,2)));
    rlength = M/rsect;
    firsttarget = zeros(M,1);
    for i = 1:rsect
        firsttarget((i-1)*rlength+1:i*rlength) = rbeamresponses((i-1)*rlength+1:i*rlength,end-i+1);
    end
    
    Azleft = struct;
    for i = 2:length(fl)
       Azleft.(fl{i}) = AzAngles.(fl{i})(AzAngles.(fl{i}) < 0);
    end
    lbeamresponses = UPA_Codebook(lambda,ElAngles(2:end),Azleft,M_V,M_H,d_V,d_H);
    lsect = 2^floor(log2(size(lbeamresponses,2)));
    llength = M/lsect;
    secondtarget = zeros(M,1);
    for i = 1:lsect
        secondtarget((i-1)*llength+1:i*llength) = lbeamresponses((i-1)*llength+1:i*llength,end-i+1);
    end
    save('Hedge32.mat','firsttarget','secondtarget');
    %% plot
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
    caxis([0 pow2db(rlength)]); % limit the color map
    colormap(flipud(hot));
    hBar = colorbar;
    set(hBar, 'TickLabelInterpreter', 'latex');
    axis equal;
    set(gca,'color',[0.9 0.9 0.9]);
    set(gca,'fontsize',18);
    view(122,30);
    hold on;
    xlim([0,1]); % To cut the sphere in half

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
    caxis([0 pow2db(llength)]); % limit the color map
    colormap(flipud(hot));
    hBar = colorbar;
    set(hBar, 'TickLabelInterpreter', 'latex');
    axis equal;
    set(gca,'color',[0.9 0.9 0.9]);
    set(gca,'fontsize',18);
    view(122,30);
    hold on;
    xlim([0,1]); % To cut the sphere in half

    
end
