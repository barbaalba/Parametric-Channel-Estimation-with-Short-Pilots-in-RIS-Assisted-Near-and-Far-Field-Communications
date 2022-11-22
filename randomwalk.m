function [x_t,y_t] = randomwalk(NumUE,RandomLength,Xmax,Ymax,SPlim)
x_t = zeros(NumUE,RandomLength);
y_t = zeros(NumUE,RandomLength);
for m=1:NumUE
  for n = 1:RandomLength % Looping all values of N into x_t(n).
      SP = unifrnd(SPlim(1),SPlim(2));
      if x_t(m,n) < (Xmax-SP) && x_t(m,n)>(-Xmax+SP)
          A = SP*sign(randn); % Generates either +1/-1 depending on the SIGN of RAND.
          x_t(m,n+1) = x_t(m,n) + A;
      else
          A = -SP * sign(x_t(m,n));
          x_t(m,n+1) = x_t(m,n) + A;
      end
      if y_t(m,n) < Ymax && y_t(m,n) > -Ymax
          A = SP*sign(randn); % Generates either +1/-1 depending on the SIGN of RAND.
          y_t(m,n+1) = y_t(m,n) + A;
      else
          A = -SP*sign(y_t(m,n));
          y_t(m,n+1) = y_t(m,n) + A;
      end
  end
end
