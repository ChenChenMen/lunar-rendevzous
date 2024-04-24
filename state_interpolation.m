function interp_result = state_interpolation(interp,queried_time)
% this function assumes provided reference trajectory bounds queried times
% also starts at the aligned first element
ref_traj = interp.ref_traj; ref_dt = interp.ref_dt;
xlinear = @(rt,k,kp1,h) rt.*(kp1-k)./h + k;

interp_result = nan(9,length(queried_time));
for i = 1:length(queried_time)
    qt = queried_time(i);
    start_ind = ceil(qt/ref_dt);
    residual_t = qt - (ceil(qt/ref_dt)-1)*ref_dt;
    if start_ind == size(ref_traj,2)
        interp_result(:,i) = ref_traj(:,end);
        break
    end
    interp_result(1:6,i) = xlinear(residual_t,ref_traj(1:6,start_ind),ref_traj(1:6,start_ind+1),ref_dt)';
    interp_result(7:9,i) = xlinear(residual_t,ref_traj(7:9,start_ind),ref_traj(7:9,start_ind+1),ref_dt)';
end
