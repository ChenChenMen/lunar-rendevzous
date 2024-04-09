function [A,B,STM] = CR2BRM(targ_orb)
% Circular Restricted 2 Body Relative Motion
% formulated in target centered LVLH frame
% x - radially outwards
% y - along the orbit track
% z - angular momentum direction

% mu, orbiting body gravitational parameter
mu = targ_orb.moon_mu;
% a, semi-major axis of the target spacecraft
a = targ_orb.a;

% n, mean motion
n = sqrt(mu/a.^3);

% state matrix
A = [zeros(3) eye(3); diag([3*n^2 0 -n^2]) [0 2*n 0; -2*n 0 0; 0 0 0]];
% input matrix
B = [zeros(3); eye(3)];

% state transition matrix
STM_rr = @(t) [4-3*cos(n*t) 0 0; 6*(sin(n*t)-n*t) 1 0; 0 0 cos(n*t)];
STM_rv = @(t) [1/n*sin(n*t) 2/n*(1-cos(n*t)) 0;
    2/n*(cos(n*t)-1) 1/n*(4*sin(n*t)-3*n*t) 0;
    0 0 1/n*sin(n*t)];
STM_vr = @(t) [3*n*sin(n*t) 0 0; 6*n*(cos(n*t)-1) 0 0; 0 0 -n*sin(n*t)];
STM_vv = @(t) [cos(n*t) 2*sin(n*t) 0; -2*sin(n*t) 4*cos(n*t)-3 0; 0 0 cos(n*t)];
STM = @(t) [STM_rr STM_rv; STM_vr STM_vv];
end