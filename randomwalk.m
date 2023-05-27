function [x_t,y_t] = randomwalk(NumUE,RandomLength,Xmax,Ymax,SPlim,d_bjo,RIS)
% We break each second into many sub period 
breakfreq = 5;
x_t = zeros(NumUE,RandomLength*breakfreq);
y_t = zeros(NumUE,RandomLength*breakfreq);
xinit = unifrnd(Xmax/4,Xmax/2,NumUE,1);% Initilaize the location of the user
yinit = unifrnd(-Ymax/2,Ymax/2,NumUE,1);% Initilaize the location of the user
x_t(:,1:breakfreq) = xinit;
y_t(:,1:breakfreq) = yinit;
for m=1:NumUE
  for n = 2:RandomLength % Looping all values of N into x_t(n).
      SP = unifrnd(SPlim(1),SPlim(2)); % The movement speed
      dir = unifrnd(-pi,pi); % The movement direction
      Xchange = SP/breakfreq*sin(dir);
      Ychange = SP/breakfreq*cos(dir);
      for l = 1:breakfreq
          Xnew = x_t(m,breakfreq*(n-1)+(l-1))+Xchange;
          Ynew = y_t(m,breakfreq*(n-1)+(l-1))+Ychange;
          newdist = sqrt((Xnew-RIS(1))^2+(Ynew-RIS(2))^2);
          if newdist > d_bjo && Xnew > 0 && Xnew < Xmax && Ynew < Ymax && Ynew > -Ymax
              x_t(m,breakfreq*(n-1)+l) = Xnew;
              y_t(m,breakfreq*(n-1)+l) = Ynew;
          else
              x_t(m,breakfreq*(n-1)+l) = x_t(m,breakfreq*(n-1)+(l-1))-Xchange;
              y_t(m,breakfreq*(n-1)+l) = y_t(m,breakfreq*(n-1)+(l-1))-Ychange;
          end
      end

  end
end

