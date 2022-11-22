function [phase] = nearFieldChanGen(x_t,y_t,z_t,U,lambda)

d_x = abs(U(1,:).' - x_t); % [M*T]
d_y = abs(U(2,:).' - y_t); % [M*T] 
d_z = abs(U(3,:).' - z_t); % [M*T]

dist = sqrt(d_x.^2+d_y.^2+d_z.^2); % distance with respect to antenna elements
phase = exp(-1i*2*pi*dist/lambda);

end