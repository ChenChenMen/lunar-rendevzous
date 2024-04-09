%% FUNCTION - Generalized RK4 solver for a dynamic equation
function record = rk4sim(dynFunc,init_cond,dt,tf,u_hist)
t0 = 0; x = init_cond; forces = zeros([9,1]);
t_hist = t0:dt:tf; x_hist = init_cond; force_hist = forces;
for i = 1:size(t_hist,2)
    t = t_hist(i); x_hist(:,i) = x; force_hist(:,i) = forces;
    if i <= size(u_hist,2) u = u_hist(:,i); else u = zeros([4,1]); end
    [x,forces] = singleRK4(dynFunc,t,dt,x,u);
end
record.t = t_hist; record.x = x_hist; record.force = force_hist;
