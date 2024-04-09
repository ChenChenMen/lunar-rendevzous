%% FUNCTION - plotter
function tcvh_plot(time,state,targ_orb,no_quiver,zoom_in)
% formulated in target centered LVLH frame (TCVH)
% origin at target spacecraft
% x - defined by cross product (roughly velocity direction)
% y - opposite of angular momentum
% z - radially inwards

% figure(Position=[0,0,1000,800]);

% % % plot the target as fixed point % % %
scatter3(0,0,0,"filled","o","MarkerFaceColor",[0 0 0]); hold on
% hab tube - 80 x 8 x 8 ft3
rad_hab = targ_orb.gateway.rad_hab;
len_hab = targ_orb.gateway.len_hab;
offset_hab = targ_orb.gateway.offset_hab;
[up_hab,bot_hab,side_hab] = xtube(rad_hab,len_hab,offset_hab);
surf(up_hab.x, up_hab.y, up_hab.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
surf(bot_hab.x, bot_hab.y, bot_hab.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
surf(side_hab.x, side_hab.y, side_hab.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);

% element tube - 120 x 8 x 8 ft3
rad_elem = targ_orb.gateway.rad_elem;
len_elem = targ_orb.gateway.len_elem;
offset_elem = targ_orb.gateway.offset_elem;
[up_elem,bot_elem,side_elem] = ytube(rad_elem,len_elem,offset_elem);
surf(up_elem.x, up_elem.y, up_elem.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
surf(bot_elem.x, bot_elem.y, bot_elem.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
surf(side_elem.x, side_elem.y, side_elem.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);

% airlock tube (assumed docking location) - 20 x 8 x 8 ft3
rad_air = targ_orb.gateway.rad_air;
len_air = targ_orb.gateway.len_air;
offset_air = targ_orb.gateway.offset_air;
[up_air,bot_air,side_air] = ytube(rad_air,len_air,offset_air);
surf(up_air.x, up_air.y, up_air.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
surf(bot_air.x, bot_air.y, bot_air.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);
surf(side_air.x, side_air.y, side_air.z, 'FaceColor',[0.5 0.5 0.5],"FaceAlpha",0.5,"EdgeColor",[0.7 0.7 0.7]);

% dock gate approach path cone
approach_len = targ_orb.gateway.trajcone_h;
approach_rad = approach_len/targ_orb.gateway.trajcone_poly(1);
contact_offset = targ_orb.gateway.trajcone_contact_rad*targ_orb.gateway.trajcone_poly(1);
approach_offset = targ_orb.current_target_us(1:3)'-targ_orb.current_target_us(4:6)'.*contact_offset;
[approach_up,approach_bot] = xcone(approach_rad,approach_len,approach_offset,1);
surf(approach_up.x, approach_up.y, approach_up.z, 'FaceColor',[0.1 0.5 0.1],"FaceAlpha",0.5,"EdgeColor",[0 0.7 0]);
surf(approach_bot.x, approach_bot.y, approach_bot.z, 'FaceColor',[0.1 0.5 0.1],"FaceAlpha",0.5,"EdgeColor",[0 0.7 0]);

% % % plot relative motion trajectory % % %
rho_si = state(:,1:3)*targ_orb.EM.scales.length*1e3;
rho_dot_si = state(:,4:6)*targ_orb.EM.scales.speed*1e3;
if ~exist("no_quiver","var") || ~no_quiver
    % plot the traveled trajectory in solid blue
    plot3(rho_si(:,1),rho_si(:,2),rho_si(:,3),"b","LineStyle","-","LineWidth",3);
    if exist("zoom_in","var") && zoom_in
        % plot only the last velocity vector
        quiver3(rho_si(end,1),rho_si(end,2),rho_si(end,3),...
            rho_dot_si(end,1),rho_dot_si(end,2),rho_dot_si(end,3),"LineWidth",1,"Color","red");
    else
        % plot the entire velocity vector trajectory
        quiver3(rho_si(:,1),rho_si(:,2),rho_si(:,3),...
            rho_dot_si(:,1),rho_dot_si(:,2),rho_dot_si(:,3),"LineWidth",1,"Color","red");
    end
else
    % plot the planned trajectory in dashed grey
    plot3(rho_si(:,1),rho_si(:,2),rho_si(:,3),"Color",[0.5 0.5 0.5],"LineStyle","--","LineWidth",1);
end
% plot features
xlabel("$x_{TCVH}$, (m)"); ylabel("$y_{TCVH}$, (m)"); zlabel("$z_{TCVH}$, (m)");
axis equal; xlimit = xlim(); ylimit = ylim(); zlimit = zlim();
% zoom in to center at the current position
if exist("zoom_in","var") && zoom_in
    zoom_size = 5;
    xlimit = rho_si(end,1)+[-zoom_size zoom_size];
    ylimit = rho_si(end,2)+[-zoom_size zoom_size];
    zlimit = rho_si(end,3)+[-zoom_size zoom_size];
end
xlim(xlimit); ylim(ylimit); zlim(zlimit);
hold off
