function interp_result = state_interpolation(interp,queried_time)
% this function assumes provided reference trajectory bounds queried times
% also starts at the aligned first element
ref_traj = interp.ref_traj; ref_dt = interp.ref_dt;
% dyn_dt = ref_dt./targ_orb.EM.scales.time;

% Hermite spline fit
% C = @(h,x) [1 0 0 0; 0 1 0 0; -3/h^2 -2/h 3/h^2 -1/h; 2/h^3 1/h^2 -2/h^3 1/h^2]*x;
% xspline = @(t,C,h) C(1,floor(t/h)+1)+C(2,floor(t/h)+1).*mod(t,h)+...
%     C(3,floor(t/h)+1).*mod(t,h).^2+C(4,floor(t/h)+1).*mod(t,h).^3;
xlinear = @(t,k,kp1,h) mod(t,h).*(kp1-k)./h + k;

interp_result = nan(9,length(queried_time));
for i = 1:length(queried_time)
    qt = queried_time(i);
    start_ind = ceil(qt/ref_dt);
    residual_t = mod(qt,ref_dt);
    % start_x = [ref_traj(1:6,start_ind)./targ_orb.EM.scales.state; ref_r];
    % end_x = [ref_traj(1:6,start_ind+1)./targ_orb.EM.scales.state; ref_r];
    % start_dxdt = ER3BRM(start_x,ref_traj(7:9,start_ind),targ_orb);
    % end_dxdt = ER3BRM(end_x,ref_traj(7:9,start_ind+1),targ_orb);
    % 
    % for state_ind = 1:6
    %     xk = start_x(state_ind); xkp1 = end_x(state_ind);
    %     xdotk = start_dxdt(state_ind); xdotkp1 = end_dxdt(state_ind);
    %     Ctraj = C(dyn_dt,[xk; xdotk; xkp1; xdotkp1]);
    %     interp_result(state_ind,i) = xspline(residual_t./targ_orb.EM.scales.time,Ctraj,dyn_dt)';
    % end
    interp_result(1:6,i) = xlinear(residual_t,ref_traj(1:6,start_ind),ref_traj(1:6,start_ind+1),ref_dt)';
    interp_result(7:9,i) = xlinear(residual_t,ref_traj(7:9,start_ind),ref_traj(7:9,start_ind+1),ref_dt)';
end
