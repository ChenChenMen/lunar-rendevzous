function [up,bot] = xcone(rad,len,offset,flip)
flip_sign = 1-2*flip;

d = -rad:2*rad/20:rad;
x_up = flip_sign*[zeros(size(d)); zeros(size(d)); len*ones(size(d)); len*ones(size(d))];
y_up = [zeros(size(d)); zeros(size(d)); d; d];
z_up = [zeros(size(d)); zeros(size(d));
    sqrt(rad^2*ones(size(d))-d.^2); zeros(size(d))];

x_bot = flip_sign*[zeros(size(d)); zeros(size(d)); len*ones(size(d)); len*ones(size(d))];
y_bot = [zeros(size(d)); zeros(size(d)); d; d];
z_bot = [zeros(size(d)); zeros(size(d));
    -sqrt(rad^2*ones(size(d))-d.^2); zeros(size(d))];

up.x = x_up+offset(1); up.y = y_up+offset(2); up.z = z_up+offset(3);
bot.x = x_bot+offset(1); bot.y = y_bot+offset(2); bot.z = z_bot+offset(3);
