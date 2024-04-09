function dxdt = ER3BRM(state,ctrl,targ_orb)
% Exact Restricted 3 Body Relative Motion
% https://arc.aiaa.org/doi/10.2514/1.A34390

% formulated in target centered LVLH frame (TCVH)
% x - defined by cross product (roughly velocity direction)
% y - opposite of angular momentum
% z - radially inwards
% defines Moon centered orbital frame (MCO)
% origin at the Moon's center of mass
% x - byrocenter pointing towards the Earth
% y - defined by cross product
% z - angular momentum

% state composite with the following
% 1:3 - rho, position vector of chaser in TCVH
% 4:6 - rhodot, velocity vector of chaser in TCVH
% 7:9 - r, position vecotor of target in MCO
% 10:12 - rdot, velocity vecotor of target in MCO
rho = state(1:3); rho_dot = state(4:6);
r = state(7:9); rdot = state(10:12);
rn = norm(r); rdotn = norm(rdot);

% control is formulated with 6 thrusters pointing at different direction
% thruster is ON/OFF command based represented by 0s and 1s
% ctrl variable should be 3x1 each valued by -1, 0, or 1
ctrl_acce = targ_orb.thrust*ctrl;

% mu, mass ratio between the two primaries
mu = targ_orb.mu;
% rem, position vector from the first primary to the second primary
% @ a circular Earth-Moon system assumption - rem is constant in i
rem = targ_orb.rem; remn = norm(rem);
% remdot, inertial time derivative rem in inertial
% @ a circular Earth-Moon system assumption - remdot is zero
remdot = targ_orb.remdot;

% omegami, MCO angular velocity wrt inertial frame
% @ a circular Earth-Moon system assumption - omegami is constant in k
omegami = targ_orb.omegami;
% omegamidot, MCO angular acceleration wrt inertial frame
% @ a circular Earth-Moon system assumption - omegamidot is zero
omegamidot = targ_orb.omegamidot;
% omegamidotdot, MCO angular jerk wrt inertial frame
% @ a circular Earth-Moon system assumption - omegamidotdot is zero
omegamidotdot = targ_orb.omegamidotdot;

% % % Target in MCO % % %
% re, position vector of target in MCO from the Earth to target
re = r+rem; ren = norm(re);
% h, angular momentum of target in MCO
h = cross(r,rdot); hn = norm(h);
% rdotdot, acceleration vector of target in MCO
rdotdot = -2*cross(omegami,rdot) - cross(omegamidot,r) - cross(omegami,cross(omegami,r)) + ...
    -mu*r/rn^3 - (1-mu)*(re/ren^3-rem/remn^3);
% rjerk, jerk vector of target in MCO
partmat1 = 1/rn^3*eye(3) - 3/rn^5*(r*r');
partmat2 = 1/ren^3*eye(3) - 3/ren^5*(re*re');
partmat3 = 1/remn^3*eye(3) - 3/remn^5*(rem*rem');
remdotmco = remdot - cross(omegami,rem);
rjerk = -2*cross(omegami,rdotdot) - cross(omegamidot,rdot) - cross(omegamidotdot,r) - cross(omegamidot,cross(omegami,r)) + ...
    cross(omegami,cross(omegamidot,r)) - cross(omegami,cross(omegami,rdot)) - mu*partmat1*rdot + ...
    -(1-mu)*(partmat2*(rdot+remdotmco)-partmat3*remdotmco);

% omegalm, TCVH angular velocity wrt MCO
omegalm = [0; -hn/rn^2; -rn/hn^2*dot(h,rdotdot)];
% omegalmdot, TCVH angular acceleration wrt MCO
rdotn = dot(r/rn,rdot); hdotn = norm(cross(r,rdotdot));
omegalmdot = [0; -1/rn*(hdotn/rn+2*rdotn*omegalm(2)); (rdotn/rn-2*hdotn/hn)*omegalm(3)-rn/hn^2*dot(h,rjerk)];

% % % Chaser in TCVH % % %
% omegali, TCVH angular velocity wrt inertial
omegali = omegalm + omegami;
% omegalidot, TCVH angular acceleration wrt inertial
omegalidot = omegalmdot + omegamidot - cross(omegalm,omegami);

% % % state update % % %
dxdt = nan(size(state));
dxdt(1:3) = rho_dot;
dxdt(4:6) = -2*cross(omegali,rho_dot) - cross(omegalidot,rho) - cross(omegali,cross(omegali,rho)) + ...
    mu*(r/rn^3-(r+rho)/norm(r+rho)^3) + (1-mu)*(re/ren^3-(re+rho)/norm(re+rho)^3) + ctrl_acce;
dxdt(7:9) = rdot;
dxdt(10:12) = rdotdot;
end