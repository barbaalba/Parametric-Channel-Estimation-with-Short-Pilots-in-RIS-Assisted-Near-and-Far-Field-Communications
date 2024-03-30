function [x_t,y_t] = randomwalk(NumUE,Simruntime,Xmax,Ymax,SPlim,d_bjo,d_NF,RIS,breakfreq)
% this function generate a user walking randomly inside a room generating
% its location within X-Y plane
% Input: 
%   RIS: location of the RIS reference element which is [0 0 0]
%   Simruntime: time in second
%   SPlim: min and max speed value
% Output: 
%   x_t,y_t: are the location of the user in time
% Room Scenario:
%       RIS_Loc = [0,0,0]
%   ----------RIS---------- x = 0
%   |                      |
%   |                      |
%   |        ROOM          | 
%   |                      |
%   |                      |
%   ----------------------- x = Xmax
% -Ymax                   Ymax  
x_t = zeros(NumUE,(Simruntime+1)*breakfreq);
y_t = zeros(NumUE,(Simruntime+1)*breakfreq);

% Initilaize the location of the user
rinit = d_NF;
anginit = unifrnd(0,pi,1);
xinit = rinit*sin(anginit); 
yinit = rinit*cos(anginit); 

% user is static for one second
x_t(:,1:breakfreq) = xinit; 
y_t(:,1:breakfreq) = yinit;

for m=1:NumUE
  for n = 2:Simruntime+1 % Channel realization instances
      SP = unifrnd(SPlim(1),SPlim(2)); % The movement speed, constant for one second
      dir = unifrnd(-pi,pi); % The movement direction, constant for one second
      Xchange = SP/breakfreq*sin(dir); % relocation in x-direction for each channel instance 
      Ychange = SP/breakfreq*cos(dir); % relocation in y-direction for each channel instance
      for l = 1:breakfreq
          Xnew = x_t(m,breakfreq*(n-1)+(l-1))+Xchange;
          Ynew = y_t(m,breakfreq*(n-1)+(l-1))+Ychange;
          newdist = sqrt((Xnew-RIS(1))^2+(Ynew-RIS(2))^2);
          % chack if the user does not exit the room or get to close to
          % RIS (not below björnson distance), if it does change the
          % moving direction
          if newdist > d_bjo && Xnew > 0 && Xnew < Xmax && Ynew < Ymax && Ynew > -Ymax
              x_t(m,breakfreq*(n-1)+l) = Xnew;
              y_t(m,breakfreq*(n-1)+l) = Ynew;
          else
              Xchange = -Xchange; Ychange = -Ychange;
              x_t(m,breakfreq*(n-1)+l) = x_t(m,breakfreq*(n-1)+(l-1))+Xchange;
              y_t(m,breakfreq*(n-1)+l) = y_t(m,breakfreq*(n-1)+(l-1))+Ychange;
          end
      end

  end
end

x_t = x_t(breakfreq+1:end);
y_t = y_t(breakfreq+1:end);

end

