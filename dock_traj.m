[log_x,log_param,fcost,fc,fceq] = dock_traj_1(state,targ_orb);
save("rendezv_traj.mat","log_x","log_param","targ_orb","fceq","fc","fcost");
function [xtraj,param,fcost,fc,fceq] = dock_traj_1(state,targ_orb)
% rendezvous using the nonlinear relative motion dynamics
% problem scales in target proximity metrics
% requires the state input to be in target proximity scale
scales.length = targ_orb.gateway.scales.length;
scales.speed = targ_orb.gateway.scales.speed;
scales.time = targ_orb.gateway.scales.time;
scales.EM.length = targ_orb.EM.scales.length;
scales.EM.speed = targ_orb.EM.scales.speed;
scales.EM.time = targ_orb.EM.scales.time;
param.scales = scales;

% trajectory optimizer setup
param.n = 200;
param.num_state = 6;
param.num_ctrl = 3;
param.num = param.num_state+param.num_ctrl;

% collect parameters for optimization use
param.scales.state = [scales.length*ones([3,1]); scales.speed*ones([3,1])];
param.scales.EM.state = [scales.EM.length*ones([3,1]); scales.EM.speed*ones([3,1])]*1e3;
param.scales.EM2Prox = param.scales.state./param.scales.EM.state;
param.scales.EM2Proxdt = (param.scales.state./scales.time)./(param.scales.EM.state./scales.EM.time);
param.x0 = state(1:6)./param.scales.EM2Prox;
param.init_state = state;
max_time = 600/scales.time;
param.init_dt = 15/scales.time;
param.target = targ_orb.current_target_us(1:3)'./(scales.length*ones([3,1]));

% construct initial optimization variable
init_x = nan(param.n*param.num+1,1);
init_x(end) = param.init_dt;
% ode propagation
er3brm_ode = @(t,x) ER3BRM(x,[0;0;0],targ_orb);
tspan = (0:init_x(end):param.n*init_x(end))*scales.time/scales.EM.time;
[~,param.zeroctrl_states] = ode45(er3brm_ode,tspan,state);
init_proxi_state = param.zeroctrl_states(2:end,1:param.num_state)./param.scales.EM2Prox';
init_proxi_ctrl = 2*rand(param.n,param.num_ctrl)-1;
init_x(1:end-1) = reshape([init_proxi_state init_proxi_ctrl],[param.num*param.n,1]);

% use the previous run results
load("rendezv_traj_poscon.mat","log_x","log_param");
interp.ref_traj = [[param.x0;0;0;0] [reshape(log_x(1:log_param.n*log_param.num_state),[log_param.n,log_param.num_state])';
    reshape(log_x(log_param.n*log_param.num_state+1:log_param.n*log_param.num),[log_param.n,log_param.num_ctrl])']];
interp.ref_dt = log_x(end);
queried_time = (1:param.n)*log_x(end)*log_param.n/param.n;
interp_result = state_interpolation(interp,queried_time);
init_x = [reshape(interp_result(1:param.num_state,:)',[param.num_state*param.n,1]);
    reshape(interp_result(param.num_state+1:param.num,:)',[param.num_ctrl*param.n,1]);
    log_x(end)*log_param.n/param.n];

% setup upper and lower bounds
low_bounds = [-1*ones(param.n*3,1); -2/scales.speed*ones(param.n*3,1); -1*ones(param.n*3,1); 0];
upp_bounds = [1*ones(param.n*3,1); 2/scales.speed*ones(param.n*3,1); 1*ones(param.n*3,1); max_time/param.n];

% setup cost and constraint functions
cost_func = @(x) cost(x,param,targ_orb);
nonc_func = @(x) nonl_constraint(x,param,targ_orb);
% initial guess cost and constraint
initcost = cost_func(init_x); [initc,initceq] = nonc_func(init_x);
% run optimization
options = optimoptions('fmincon','Display','iter-detailed','Algorithm','interior-point',...
    "SubproblemAlgorithm","cg",'MaxFunctionEvaluations',1e8,'MaxIterations',1000);
xtraj = fmincon(cost_func,init_x,[],[],[],[],low_bounds,upp_bounds,nonc_func,options);
% xtraj(param.n*param.num_state+1:param.n*param.num) = ctrl_cond(xtraj,param);
fcost = cost_func(xtraj); [fc,fceq] = nonc_func(xtraj);

    function val = cost(x,param,targ_orb)
        n = param.n; num = param.num; num_state = param.num_state;
        % design weighting
        w_lat_dist = 10; w_apr_dist = 1; fw_dist = 100;
        w_lat_velo = 10; w_apr_velo = 10; fw_velo = 1000;

        states = reshape(x(1:n*num_state),[n,num_state]);
        controls = ctrl_cond(x,param); dt = x(end);
        dist = (states(:,1:3)'-param.target).*param.scales.state(1:3);
        velo = (states(:,4:6)'-zeros(3,1)).*param.scales.state(4:6);
        % cost function design
        val = 20*sum(abs(controls))*dt+...
            fw_dist*norm(dist(:,end)) + fw_velo*norm(velo(:,end));
    end

    function [c,ceq] = nonl_constraint(x,param,targ_orb)
        persistent func_count;
        if isempty(func_count) 
            func_count = 0; 
        end

        % dynamics feasibility by collocation method
        n = param.n; num = param.num;
        num_state = param.num_state; num_ctrl = param.num_ctrl;
        x0 = param.x0; init_state = param.init_state; dt = x(end);
        states = [x0 reshape(x(1:num_state*n),[n,num_state])'];
        % condition the control input
        contrl_signal = x(param.num_state*param.n+1:param.n*param.num);
        contrl = [zeros(num_ctrl,1) reshape(contrl_signal,[n,num_ctrl])'];

        dxdt = nan(6,n+1);
        for i = 1:n+1
            orb_state = param.zeroctrl_states(floor(n/2),7:12)';
            if (i-1)*x(end) <= n*param.init_dt
                bot_index = min(n,floor((i-1)*x(end)/param.init_dt));
                if bot_index == 0
                    bot_orb_state = init_state(7:12);
                else
                    bot_orb_state = param.zeroctrl_states(bot_index,7:12)';
                end
                top_index = min(n,max(1,ceil((i-1)*x(end)/param.init_dt)));
                % orbital state interpolation
                orb_state = bot_orb_state+mod((i-1)*x(end),param.init_dt)*...
                    (param.zeroctrl_states(top_index,7:12)'-bot_orb_state);
            end
            dyn_states = [states(:,i).*param.scales.EM2Prox; orb_state];
            result = ER3BRM(dyn_states,contrl(:,i),targ_orb);
            dxdt(:,i) = result(1:6)./param.scales.EM2Proxdt;
        end

        xcol = 1/2*(states(:,1:end-1)+states(:,2:end))-dt/8*diff(dxdt,1,2);
        xcoldot = 3/2/dt*diff(states,1,2)-1/4*(dxdt(:,1:end-1)+dxdt(:,2:end));
        ucol = 1/2*(contrl(:,1:end-1)+contrl(:,2:end));

        dxdt_col = nan(6,n);
        for i = 1:n
            orb_state = param.zeroctrl_states(floor(n/2),7:12)';
            if (i-0.5)*x(end) <= n*param.init_dt
                bot_index = floor((i-0.5)*x(end)/param.init_dt);
                if bot_index == 0
                    bot_orb_state = init_state(7:12);
                else
                    bot_orb_state = param.zeroctrl_states(bot_index,7:12)';
                end
                top_index = max(ceil((i-0.5)*x(end)/param.init_dt));
                % orbital state interpolation
                orb_state = bot_orb_state+mod((i-0.5)*x(end),param.init_dt)*...
                    (param.zeroctrl_states(top_index,7:12)'-bot_orb_state);
            end
            dyn_states = [xcol(:,i).*param.scales.EM2Prox; orb_state];
            result = ER3BRM(dyn_states,ucol(:,i),targ_orb);
            dxdt_col(:,i) = result(1:6)./param.scales.EM2Proxdt;
        end
        ceq = reshape(xcoldot-dxdt_col,[6*n,1]);
        % end state constraint reaching target
        ceq = [ceq; states(:,end)-[param.target; zeros(3,1)]];
        % ceq = [ceq; x(num_state*n+1:n*num).*(x(num_state*n+1:n*num)-1).*(x(num_state*n+1:n*num)+1)];

        c = targ1_constraint(states,param,targ_orb);
        % c = [c; vecnorm(states(4:6,:))'-2/param.scales.time];

        % visualization
        if mod(func_count,5000) == 0
            plot_state = (states.*param.scales.EM2Prox)';
            targ_orb.axes(1) = subplot(2,2,1);
            tcvh_plot([],plot_state,targ_orb);
            rho_si = plot_state(:,1:3)*targ_orb.EM.scales.length*1e3;
            rho_dot_si = plot_state(:,4:6)*targ_orb.EM.scales.speed*1e3;
            targ_orb.axes(2) = subplot(2,2,2);
            time = dt*(1:size(states,2))*param.scales.time;
            plot(time,rho_si(:,1),"k","LineWidth",2); hold on
            plot(time,rho_si(:,2),"b","LineWidth",2); grid on
            plot(time,rho_si(:,3),"r","LineWidth",2); hold off
            legend("$\rho_x$","$\rho_y$","$\rho_z$","interpreter","latex");
            xlabel("Time, s"); ylabel("rel. distance, m");

            targ_orb.axes(2) = subplot(2,2,3);
            plot(time,rho_dot_si(:,1),"k--","LineWidth",2); hold on
            plot(time,rho_dot_si(:,2),"b--","LineWidth",2); grid on
            plot(time,rho_dot_si(:,3),"r--","LineWidth",2); hold off
            legend("$\delta v_x$","$\delta v_y$","$\delta v_z$","interpreter","latex");
            xlabel("Time, s"); ylabel("rel. velocity, m/s");

            targ_orb.axes(4) = subplot(2,2,4);
            plot(time(1:end-1),xcoldot'-dxdt_col'); grid on; hold off
            xlabel("Time, s"); ylabel("collocation error");
            drawnow
        end
        func_count = func_count + 1;
    end

    function real_ctrl = ctrl_cond(x,param)
        real_ctrl = x(param.num_state*param.n+1:param.n*param.num);
        zero_ind = abs(real_ctrl)<1e-6; real_ctrl(zero_ind) = 0;
        real_ctrl(~zero_ind) = real_ctrl(~zero_ind)./abs(real_ctrl(~zero_ind));
    end

end
