function c = targ1_constraint(state,param,targ_orb)
% Docking target 1 constraint - Orion Capsule
% input state unit should be proximity scaled
target_pos = targ_orb.current_target_us(1:3)'./param.scales.state(1:3);

% approach position cone
trajcone_poly = targ_orb.gateway.trajcone_poly;
trajcone_h = targ_orb.gateway.trajcone_h/param.scales.length;
traj_contact_offset = targ_orb.gateway.trajcone_contact_rad*targ_orb.gateway.trajcone_poly(1);
rho_dock = state(1:3,:) - (target_pos - targ_orb.current_target_us(4:6)'.*traj_contact_offset/param.scales.length);
% position cone constraint
r = vecnorm(rho_dock(2:3,:)); % y,z direction
h = rho_dock(1,:); h_ref = -h;
trajcone_order = length(trajcone_poly)-1;
trajcone = trajcone_poly(end) + zeros(size(r));
r_curr_order = r;
for i = 1:trajcone_order
    trajcone = trajcone + trajcone_poly(i).*r_curr_order;
    r_curr_order = r_curr_order.*r;
end
below_ind = h_ref<trajcone_h;
c_h = [trajcone_h-h_ref(~below_ind) trajcone(below_ind)-h_ref(below_ind)]';

% approach x velocity steps
velostep_h = targ_orb.gateway.velostep/param.scales.length;
out_ind = h_ref>=velostep_h(1);
v_constraint = repmat(velostep_h(1)-h_ref(out_ind),[3,1]);
for i = 1:length(velostep_h)-1
    ind = h_ref<velostep_h(i) & h_ref>=velostep_h(i+1);
    ind_size = sum(ind);
    limit_v = targ_orb.gateway.velolimit(:,i:i+1)/param.scales.speed;
    v_constraint = [v_constraint abs(state(4:6,ind))-(limit_v(:,1)-(limit_v(:,1)-limit_v(:,2))*(velostep_h(i)-h_ref(ind))/(velostep_h(i)-velostep_h(i+1)))];
end
c_v = reshape(v_constraint',[3*(param.n+1),1]);
c = [c_h; c_v];
end
