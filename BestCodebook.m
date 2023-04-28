clc;clear;close all;
%% Env Initialization
freq = 28e9; % Central frequency
lambda = physconst('LightSpeed') / freq; % Wavelength
M_H = 8; M_V = 8; M = M_H*M_V;
d_H = 1/2; d_V = 1/2; %In wavelengths
Th = -pi/2:4*pi/180:pi/2;
Az = -pi/2:4*pi/180:pi/2;
Num = zeros(length(Th),length(Az));
for thidx = 1:length(Th)
    parfor azidx = 1:length(Az)
        [ElAngles,AzAngles,CBL] = UPA_BasisElupnew(M_V,M_H,d_V,d_H,Az(azidx),Th(thidx));
        beamresponses = UPA_Codebook(lambda,ElAngles,AzAngles,M_V,M_H,d_V,d_H);
        x = abs(1/M*beamresponses'*beamresponses);
        orthnum = sum(x>0.5,'all');
%         disp(['Number of orthogonal basis: ' num2str(orthnum)]);
%         if CBL~= orthnum && CBL+2~=orthnum
%             disp(['az = ', num2str(Az(azidx)), ' El = ', num2str(Th(thidx))]);
%             pause(5);
%         end
        if CBL~= orthnum
            diff = abs(orthnum - CBL);
            CBL = CBL - diff - diff/2;
        end
        Num(thidx,azidx) = CBL;
    end
end
[X,Y] = meshgrid(-pi/2:4*pi/180:pi/2,-pi/2:4*pi/180:pi/2);
s = [-45 20];
sl = surfl(X,Y,Num,s,'light');
sl(2).Style = 'ambient';
