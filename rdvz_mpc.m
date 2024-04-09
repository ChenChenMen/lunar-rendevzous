function ctrl = rdvz_mpc(state,ref_traj,targ_orb)
% rendezvous using the nonlinear relative motion dynamics
% controller setup
if ~isfield(targ_orb,"mpc_horizon_steps")
    param.n = size(ref_traj,2);
else
    param.n = targ_orb.mpc_horizon_steps;
end
param.ref_traj = ref_traj;
param.n_state = 6; param.n_ctrl = 3;
param.ntot = param.n_state+param.n_ctrl;
param.dt = targ_orb.dt;
param.ref_r = state(7:12);

% real SI TCVH state
unscaled_rho = state(1:6).*targ_orb.EM.scales.state;

% % extract segments of reference trajectory
% [~,start_ind] = min(vecnorm(targ_orb.ref_traj(1:3,:)-unscaled_rho(1:3)));
% range_ind = ceil(param.dt*param.n/targ_orb.ref_dt);
% exceed_ind = start_ind + range_ind - size(targ_orb.ref_traj,2);
% if exceed_ind > 0
%     raw_ref_traj = [targ_orb.ref_traj(:,start_ind:end) repmat(targ_orb.ref_traj(:,end),[1,exceed_ind])];
% else
%     raw_ref_traj = targ_orb.ref_traj(:,start_ind:start_ind+range_ind);
% end
% % interpolation - assume the start time is aligned
% interp.ref_traj = raw_ref_traj; interp.ref_dt = targ_orb.ref_dt;
% param.ref_traj = state_interpolation(interp,param.ref_r,(0:(param.n-1))*param.dt,targ_orb);

% construct initial optimization variable
delta_x = repmat(unscaled_rho,[1,param.n])-param.ref_traj(1:param.n_state,:);
delta_u = zeros([param.n_ctrl,param.n]); %-param.ref_traj(param.n_state+1:param.ntot,:);
init_x = [reshape(delta_x,[param.n*param.n_state,1]); reshape(delta_u,[param.n*param.n_ctrl,1])];

% linear state propagation equality constraint
[A,B] = ERL3BRM(param.ref_r,targ_orb);
comb_mat = [A B; zeros(param.n_ctrl,param.ntot)];
comb_mat_d = expm(comb_mat*param.dt./targ_orb.EM.scales.time);
Ad = diag(targ_orb.EM.scales.state)*comb_mat_d(1:param.n_state,1:param.n_state)*diag(1./targ_orb.EM.scales.state);
Bd = diag(targ_orb.EM.scales.state)*comb_mat_d(1:param.n_state,param.n_state+1:end);

sparse_a_mat = [zeros(1,param.n); [eye(param.n-1) zeros(param.n-1,1)]];
Aeq = [-eye(param.n*param.n_state)+kron(sparse_a_mat,Ad) kron(eye(param.n),Bd)];
Beq = [-Ad*(unscaled_rho-param.ref_traj(1:param.n_state,1)); zeros((param.n-1)*param.n_state,1)];

% setup upper and lower bounds
low_bounds = [-Inf*ones(param.n*param.n_state,1); -1*ones(param.n*param.n_ctrl,1)];
upp_bounds = [Inf*ones(param.n*param.n_state,1); 1*ones(param.n*param.n_ctrl,1)];

% setup cost and constraint functions
cost_func = @(x) cost(x,param,targ_orb);
% nonl_func = @(x) nonl_constraint(x,param);

% initial condition cost
init_cost = cost_func(init_x);
disp("Initial Cost is: "+num2str(init_cost));

options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1000,'EnableFeasibilityMode',true);
x = fmincon(cost_func,init_x,[],[],Aeq,Beq,low_bounds,upp_bounds,[],options);
sol_contrl = reshape(x(param.n*param.n_state+1:param.n*param.ntot),[param.n_ctrl,param.n])+param.ref_traj(param.n_state+1:param.ntot,:);
ctrl = sol_contrl(1:3,1);

    function val = cost(x,param,targ_orb)
        persistent func_count;
        if isempty(func_count)
            func_count = 0; 
        end

        states = reshape(x(1:param.n*param.n_state),[param.n_state,param.n]);
        contrl = reshape(x(param.n*param.n_state+1:param.n*param.ntot),[param.n_ctrl,param.n])';
        val = sum(abs(contrl),"all") + 10*sum(vecnorm(states(1:3,:))) + sum(vecnorm(states(4:6,:)));

        % visualization
        % if mod(func_count,5000) == 0
        %     plot_state = states + param.ref_traj(1:6,:);
        %     % targ_orb.axes(2) = subplot(2,2,2);
        %     % plot3(param.ref_traj(1,:),param.ref_traj(2,:),param.ref_traj(3,:),"k","LineWidth",2); hold on
        %     % plot3(plot_state(1,:),plot_state(2,:),plot_state(3,:),"b","LineWidth",2); grid on; hold off;
        %     targ_orb.axes(3) = subplot(2,2,3);
        %     plot3(param.ref_traj(1,:),param.ref_traj(2,:),param.ref_traj(3,:),"k","LineWidth",2); hold on
        %     plot3(plot_state(1,:),plot_state(2,:),plot_state(3,:),"b","LineWidth",2); grid on; hold off;
        %     targ_orb.axes(4) = subplot(2,2,4);
        %     plot3(param.ref_traj(4,:),param.ref_traj(5,:),param.ref_traj(6,:),"k","LineWidth",2); hold on
        %     plot3(plot_state(4,:),plot_state(5,:),plot_state(6,:),"b","LineWidth",2); grid on; hold off;
        %     drawnow
        % end
        func_count = func_count + 1;
    end

    % function [c,ceq] = nonl_constraint(x,param)
    %     contrl = reshape(x(param.n*param.n_state+1:param.n*param.ntot),[param.n_ctrl,param.n])';
    %     real_ctrl = reshape(contrl+param.ref_traj(param.n_state+1:param.ntot,:)',[param.n_ctrl*param.n,1]);
    %     c = [real_ctrl.*(0.99-real_ctrl); real_ctrl.*(0.99+real_ctrl);]; ceq = [];
    % end

end
