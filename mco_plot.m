%% FUNCTION - plotter
function mco_plot(time, state, targ_orb)
% origin at the Moon's center of mass
% x - byrocenter pointing towards the Earth
% y - defined by cross product
% z - angular momentum

R_moon = targ_orb.radi_moon/targ_orb.EM.scales.length;
R_earth = targ_orb.radi_earth/targ_orb.EM.scales.length;

figure(Position=[0,0,1000,800]);
% % % plot the Moon % % %
[X,Y,Z] = sphere(100);
X = X*R_moon; Y = Y*R_moon; Z = Z*R_moon;
surf(X,Y,Z,"FaceColor",[0.5 0.5 0.5],"EdgeColor","none"); hold on

% % % % plot the Earth % % %
[X,Y,Z] = sphere(100);
X = X*R_earth + norm(targ_orb.rem); Y = Y*R_earth; Z = Z*R_earth;
surf(X,Y,Z,"FaceColor",[115 215 255]/255,"EdgeColor","none"); grid on

% % % plot the target trajectory % % %
plot3(state(:,7),state(:,8),state(:,9),"LineWidth",2,"Color",[0.4 0.4 0.4]);

% % % plot the chaser trajectory % % %
chaser_MCO = state(:,7:9) + state(:,1:3);
plot3(chaser_MCO(:,1),chaser_MCO(:,2),chaser_MCO(:,3),"r","LineWidth",2);

% % % figure parameters % % %
xlabel("x, DU"); ylabel("y, DU"); zlabel("z, DU");
legend("Moon","Earth","Chief S/C","Deputy S/C");
axis equal;
