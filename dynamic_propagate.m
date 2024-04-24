clear; clc; close all
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor',[1,1,1])
set(groot,'defaultAxesFontSize',12)

do_ctrl = false;

%% provide system specs
targ_orb.dt = 0.1;
targ_orb.mu = 0.01213;
% Earth Moon system setup
targ_orb.rem = [-1; 0; 0];
targ_orb.remdot = zeros([3,1]);
targ_orb.omegami = [0; 0; 1];
targ_orb.omegamidot = zeros([3,1]);
targ_orb.omegamidotdot = zeros([3,1]);
% Scales setup
targ_orb.EM.scales.time = 4.348*24*3600;
targ_orb.EM.scales.length = 384400; % km
targ_orb.radi_moon = 1740; % km
targ_orb.radi_earth = 6378; % km
% thruster specs - pulse thruster 1m/s^2 specific force
targ_orb.EM.scales.acce = targ_orb.EM.scales.length/targ_orb.EM.scales.time^2;
targ_orb.thrust = 1/1e3/targ_orb.EM.scales.acce;
% specified docking location unscaled
targ_orb.target_list = [
    -70*0.3048 0 0 -1 0 0;
    -30*0.3048 -20*0.3048 0 0 -1 0];
targ_orb.current_target_us = targ_orb.target_list(1,:);
% space station abstracted dimensions
rad_hab = 8*0.3048; len_hab = 100*0.3048; offset_hab = [-70*0.3048; 0; 0];
targ_orb.gateway.rad_hab = rad_hab;
targ_orb.gateway.len_hab = len_hab;
targ_orb.gateway.offset_hab = offset_hab;
rad_elem = 8*0.3048; len_elem = 120*0.3048; offset_elem = [0; -len_elem/2; 0];
targ_orb.gateway.rad_elem = rad_elem;
targ_orb.gateway.len_elem = len_elem;
targ_orb.gateway.offset_elem = offset_elem;
rad_air = 8*0.3048; len_air = 20*0.3048; offset_air = [-30*0.3048; -len_air; 0];
targ_orb.gateway.rad_air = rad_air;
targ_orb.gateway.len_air = len_air;
targ_orb.gateway.offset_air = offset_air;
% proximity safety parameters
targ_orb.gateway.collision_free_zone_thickness = 5; % in m
targ_orb.gateway.trajcone_contact_rad = 1*0.3048;
targ_orb.gateway.trajcone_h = 80; % in m
targ_orb.gateway.trajcone_poly = [10 0];
targ_orb.gateway.velostep = [80 40 20 10 5 0]; 
targ_orb.gateway.velolimit = [
    2 1 0.5 0.25 0.1 0.1;
    0.5 0.25 0.1 0.05 0.02 0.02;
    0.5 0.25 0.1 0.05 0.02 0.02];
% deputy spacecraft abstracted shape
a = -pi:pi/2:pi; ph = pi/4;
dpt.side.x = [cos(a+ph); cos(a+ph)]/cos(ph);
dpt.side.y = [sin(a+ph); sin(a+ph)]/sin(ph);
dpt.side.z = [-ones(size(a)); ones(size(a))];
dpt.top.x = [1 -1; 1 -1]; dpt.top.y = [1 1; -1 -1]; dpt.top.z = [1 1; 1 1];
dpt.bot.x = [1 -1; 1 -1]; dpt.bot.y = [1 1; -1 -1]; dpt.bot.z = [-1 -1; -1 -1];
% deputy spacecraft thruster location
dpt.loc_thrust = [-1*eye(3) eye(3)];
dpt.point_thrust = [-1*eye(3) -1*eye(3)];

%% initial condition setup
targ_orb.apo_alti = 3000; % km
targ_orb.apo_velo = -1.408; % km/s
targ_orb.EM.scales.speed = targ_orb.EM.scales.length/targ_orb.EM.scales.time;
% chaser in TCVH
rho = [-250; 20; 30]/1e3; % m (km)
rho_dot = [0; 0; 0]/1e3; % m/s (km/s)
chaser_rTCVH = rho/targ_orb.EM.scales.length;
chaser_rdotTCVH = rho_dot/targ_orb.EM.scales.speed;
targ_orb.gateway.scales.length = norm(rho)*1e3;
targ_orb.gateway.scales.speed = norm(rho)*1e3;
targ_orb.gateway.scales.time = 1;
if norm(rho_dot) == 0
    targ_orb.gateway.scales.time = targ_orb.gateway.scales.length/rad_hab;
    targ_orb.gateway.scales.speed = targ_orb.gateway.scales.length/targ_orb.gateway.scales.time;
end
% target in MCO
rz = (targ_orb.apo_alti+targ_orb.radi_moon)/targ_orb.EM.scales.length;
target_rMCO = [0; 0; rz];
target_rdotMCO = [0; targ_orb.apo_velo; 0]/targ_orb.EM.scales.speed;
state = [chaser_rTCVH; chaser_rdotTCVH; target_rMCO; target_rdotMCO];
targ_orb.EM.scales.state = 1e3*[targ_orb.EM.scales.length*ones(3,1); targ_orb.EM.scales.speed*ones(3,1)];

%% load trajectory file if exist
traj_file_name = "rendezv_traj_poscon.mat";
if exist(traj_file_name,"file")
    load(traj_file_name,"log_param","log_x");
    log_traj = reshape(log_x(1:log_param.n*log_param.num),[log_param.n,log_param.num]);
    targ_orb.ref_traj = [log_traj(:,1:6)'.*log_param.scales.state; log_traj(:,7:9)'];
    targ_orb.ref_dt = log_x(end).*log_param.scales.time;
    do_ctrl = true;
else
    % standard ODE propagation
    tpropg = [0 5*3600]/targ_orb.EM.scales.time; % in TU
    er3brm_ode = @(t,x) ER3BRM(x,[0;0;0],targ_orb);
    [t,sol_state] = ode89(er3brm_ode,tpropg,state);
end

%% interpolate the loaded trajectory
targ_orb.mpc_horizon_steps = 10;
if do_ctrl
    fig = figure('position', [0, 0, 1920, 1200]);
    total_time = targ_orb.ref_dt*log_param.n;
    t = targ_orb.dt:targ_orb.dt:total_time;
    init_state = [rho*1e3; rho_dot*1e3; zeros(3,1)];
    interp.ref_traj = [init_state targ_orb.ref_traj]; interp.ref_dt = targ_orb.ref_dt;
    ref_traj = [state_interpolation(interp,t) repmat(targ_orb.ref_traj(:,end),[1,targ_orb.mpc_horizon_steps])];

    %% controller deployment with reference trajectory
    % % customized ODE propagation with controls
    dynFunc = @(t,x,u) ER3BRM(x,u,targ_orb);
    for i = 1:length(t)
        range_ind = i+targ_orb.mpc_horizon_steps-1-size(ref_traj,2);
        time = t(i); future_ref_traj = ref_traj(:,i:i+targ_orb.mpc_horizon_steps-1);
        % compute the control input from MPC with actuation noise
        ctrl = rdvz_mpc(state,future_ref_traj,targ_orb).*(1+0.03*rand([3,1]));
        ctrl = max(-1, min(ctrl,1));

        % record the states and controls
        us_state = state(1:6).*targ_orb.EM.scales.state;
        t(i) = time; sol_us_state(:,i) = us_state;
        sol_state(i,:) = state'; sol_ctrl(i,:) = ctrl';

        % plot the current propagation with reference trajectory
        targ_orb.axes(1) = subplot(2,3,2);
        tcvh_plot(t,(targ_orb.ref_traj(1:6,:)./targ_orb.EM.scales.state)',targ_orb,true);
        hold on; tcvh_plot(t,sol_state,targ_orb,false,true);
        title("Deputy S/C centered view"); hold off

        targ_orb.axes(2) = subplot(2,3,1);
        tcvh_plot(t,(targ_orb.ref_traj(1:6,:)./targ_orb.EM.scales.state)',targ_orb,true);
        hold on; tcvh_plot(t,sol_state,targ_orb,false,false);
        title("Chief S/C centered view"); hold off

        targ_orb.axes(3) = subplot(2,3,4);
        plot(t,ref_traj(1,1:length(t)),"k-","LineWidth",1); hold on
        plot(t,ref_traj(2,1:length(t)),"k--","LineWidth",1); grid on
        plot(t,ref_traj(3,1:length(t)),"k-.","LineWidth",1);
        plot(t(1:i),sol_us_state(1,:),"r-","LineWidth",3);
        plot(t(1:i),sol_us_state(2,:),"g--","LineWidth",3);
        plot(t(1:i),sol_us_state(3,:),"b-.","LineWidth",3);
        xlabel("t"); ylabel("x, m");
        legend("","","","x TCVH","y TCVH","z TCVH","interpreter","latex","Location","southeast");
        title("Time history of position tracking"); hold off

        targ_orb.axes(4) = subplot(2,3,5);
        state_diff = sol_us_state(1:3,:)-ref_traj(1:3,1:i);
        plot(t(1:i),state_diff(1,:),"r-","LineWidth",1); hold on
        plot(t(1:i),state_diff(2,:),"g--","LineWidth",1); grid on
        plot(t(1:i),state_diff(3,:),"b-.","LineWidth",1);
        xlabel("t"); ylabel("$\delta x$, m");
        limit_y = ylim(); ylim([min(-0.15,limit_y(1)) max(0.15,limit_y(2))]);
        legend("x TCVH","y TCVH","z TCVH","interpreter","latex","Location","southeast");
        title("Time history of tracking error"); hold off

        targ_orb.axes(5) = subplot(2,3,3);
        view_range = 4;
        real_ctrl = dpt.point_thrust.*[max(0,ctrl); min(0,ctrl)]';
        surf(dpt.top.x,dpt.top.y,dpt.top.z,"FaceColor",[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]); hold on
        surf(dpt.bot.x,dpt.bot.y,dpt.bot.z,"FaceColor",[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]); grid on
        surf(dpt.side.x,dpt.side.y,dpt.side.z,"FaceColor",[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
        quiver3(dpt.loc_thrust(1,:),dpt.loc_thrust(2,:),dpt.loc_thrust(3,:),...
            real_ctrl(1,:),real_ctrl(2,:),real_ctrl(3,:),"r","LineWidth",3);
        axis equal; xlabel("$x_{TCVH}$, (m)"); ylabel("$y_{TCVH}$, (m)"); zlabel("$z_{TCVH}$, (m)"); view(-45,15)
        xlim([-view_range view_range]); ylim([-view_range view_range]); zlim([-view_range view_range]);
        title("Deputy S/C thruster actuation"); hold off

        targ_orb.axes(6) = subplot(2,3,6);
        plot(t(1:i),sol_ctrl(1:i,1)*100,"r-","LineWidth",1); hold on
        plot(t(1:i),sol_ctrl(1:i,2)*100,"g--","LineWidth",1); grid on
        plot(t(1:i),sol_ctrl(1:i,3)*100,"b-.","LineWidth",1);
        xlabel("t"); ylabel("u, $\%$"); ylim([-100 100]);
        legend("x TCVH","y TCVH","z TCVH","interpreter","latex","Location","northeast");
        title("Time history of control actuation"); hold off

        drawnow; F(i) = getframe(fig);

        % exact dynamics propagation
        state = singleRK4(dynFunc,time,targ_orb.dt/targ_orb.EM.scales.time,state,ctrl);
    end

    vidname = "mpc_propagate.mp4";
    v = VideoWriter(vidname,"MPEG-4"); open(v);
    for k = 1:length(F)
        writeVideo(v,F(k)); end
    close(v);

    %% data analysis and record
    save("MPC_propagation.mat","t","sol_us_state","sol_state","sol_ctrl","traj_file_name");
end

mco_plot(t,sol_state,targ_orb);
figure(Position=[0,0,1000,800]);
tcvh_plot(t,sol_state,targ_orb);