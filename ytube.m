function [up,bot,side] = ytube(rad,len,offset)
d = -rad:2*rad/20:rad;
x_up = [d; d; d; d];
y_up = [zeros(size(d)); zeros(size(d)); len*ones(size(d)); len*ones(size(d))];
z_up = [zeros(size(d)); sqrt(rad^2*ones(size(d))-d.^2);
    sqrt(rad^2*ones(size(d))-d.^2); zeros(size(d))];

x_bot = [d; d; d; d];
y_bot = [zeros(size(d)); zeros(size(d)); len*ones(size(d)); len*ones(size(d))];
z_bot = [-sqrt(rad^2*ones(size(d))-d.^2); zeros(size(d));
    zeros(size(d)); -sqrt(rad^2*ones(size(d))-d.^2)];

x_side = [d; d];
y_side = [zeros(size(d)); len*ones(size(d))];
z_side = [-sqrt(rad^2*ones(size(d))-d.^2); -sqrt(rad^2*ones(size(d))-d.^2)];

up.x = x_up+offset(1); up.y = y_up+offset(2); up.z = z_up+offset(3);
bot.x = x_bot+offset(1); bot.y = y_bot+offset(2); bot.z = z_bot+offset(3);
side.x = x_side+offset(1); side.y = y_side+offset(2); side.z = z_side+offset(3);
