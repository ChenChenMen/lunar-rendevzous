function c = targ1_constraint(state,param,targ_orb)
% Docking target 1 constraint - Orion Capsule
% input state unit should be proximity scaled
target_pos = targ_orb.current_target_us(1:3)'./param.scales.state(1:3);
cone_poly = targ_orb.gateway.trajcone_poly;
cone_h = targ_orb.gateway.trajcone_h/param.scales.length;
contact_offset = targ_orb.gateway.trajcone_contact_rad*targ_orb.gateway.trajcone_poly(1);
rho_dock = state(1:3,:) - (target_pos - targ_orb.current_target_us(4:6)'.*contact_offset/param.scales.length);
% cone constraint
r = vecnorm(rho_dock(2:3,:)); % y,z direction
h = rho_dock(1,:); h_ref = -h;
cone_order = length(cone_poly)-1;
cone = cone_poly(end) + zeros(size(r));
r_curr_order = r;
for i = 1:cone_order
    cone = cone + cone_poly(i).*r_curr_order;
    r_curr_order = r_curr_order.*r;
end
below_ind = h_ref<cone_h;
c = [cone_h-h_ref(~below_ind) cone(below_ind)-h_ref(below_ind)]';
end
