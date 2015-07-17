options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;
options.use_bullet = false;
options.enable_fastqp = false;
options.ignore_self_collisions = true;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('Planar_CBSE_Window.URDF');

G = -10;
p = PlanarRigidBodyManipulator(urdf, options);
%p = RigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt, options);

tf = 1;
x0 = [0; 0.12; pi/3; 0; 0; 0];
% x0 = [0; 0; 1; pi/6; 0; 0; 0; 0; 0; 0; 0; 0];

xtraj_ts = simulate(r, [0 tf], x0);

v = r.constructVisualizer;
% v = v.enableLCMGLInertiaEllipsoids();
v.playback(xtraj_ts, struct('slider', true));

traj = xtraj_ts.eval(xtraj_ts.tt);
x = traj(1, :);
z = traj(2, :);
theta = traj(3, :);
xdot = traj(4, :);
zdot = traj(5, :);
thetadot = traj(6, :);

% x = traj(1, :);
% y = traj(2, :);
% z = traj(3, :);
% roll = traj(4, :);
% xdot = traj(7, :);
% zdot = traj(9, :);
% rolldot = traj(10, :);
% pitchdot = traj(11, :);
% yawdot = traj(12, :);

times = linspace(0, tf, length(x));
dt = diff(times);
sensor_inds = round(linspace(1, numel(xtraj_ts.tt), 120*tf));
sensor_inds = sensor_inds(1:50);
tempxdot = xdot(sensor_inds);
tempzdot = zdot(sensor_inds);
temptimes = times(sensor_inds);
xddot = diff(tempxdot)./diff(temptimes);
zddot = diff(tempzdot)./diff(temptimes);
theta_dot = thetadot(sensor_inds);
sensor_inds = sensor_inds(1:numel(sensor_inds)-1);

for i = 1:numel(xddot)
   j = sensor_inds(i);
   ang = theta(j);
   R = [cos(ang), sin(ang); -sin(ang), cos(ang)];
   val = R*[xddot(i); zddot(i)];
   xddot(i) = val(1);
   zddot(i) = val(2);
end

N = numel(sensor_inds);

times = times(1:numel(times)-1) - times(1);
times = times';

x0 = awgn(x0, 25);
% x0 = rand(6, 1);

noisy_thetadot = theta_dot(2:end);

%add noise model
% [b, a] = butter(randi(10), rand);
% xddot = filtfilt(b, a, xddot);
% zddot = filtfilt(b, a, zddot);
% noisy_thetadot = filtfilt(b, a, noisy_thetadot);

inds = round(linspace(1, size(times, 1), N));

v = r.constructVisualizer;
v.display_dt = 0.01;
poses = zeros(r.getNumStates(), numel(x));
poses([1, 2, 3], :) = [x; z; theta];
% poses([1, 3, 4], :) = [x; z; roll];
dttimes = linspace(times(1), times(end), numel(x));
xtraj_constructed = DTTrajectory(dttimes, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
% v = v.enableLCMGLInertiaEllipsoids();
% v.draw(0, x0);
v.playback(xtraj_constructed, struct('slider', true));

% figure
% plot(times(inds), theta(inds), '*');
traj_sim = xtraj_ts;
ts_sim = traj_sim.getBreaks();
traj_sim = PPTrajectory(foh(ts_sim,traj_sim.eval(ts_sim)));
ts = times(sensor_inds);
non_meas = zeros(numel(xddot), 1);
data = [ts, xddot', non_meas, zddot', non_meas, noisy_thetadot', non_meas];
[r, xtraj, utraj, ltraj, zz, F, info, prog] = contactBasedStateEstimator(r, N, x0, data, traj_sim);

traj_ts = xtraj.getBreaks();
traj_eval = xtraj.eval(traj_ts);
x_calc = traj_eval(1, :);
z_calc = traj_eval(2, :);
theta_calc = traj_eval(3, :);
xdot_calc = traj_eval(4, :);
zdot_calc = traj_eval(5, :);
thetadot_calc = traj_eval(6, :);

% onlz plot the results of integration here
plot_inds = round(linspace(1, numel(xtraj_ts.tt), 120*tf));
plot_inds = plot_inds(2:50);
plotx = x(plot_inds);
plotz = z(plot_inds);
plottheta = theta(plot_inds);
plotxdot = xdot(plot_inds);
plotzdot = zdot(plot_inds);
plotthetadot = thetadot(plot_inds);

figure
subplot(3, 2, 1)
plot(times(sensor_inds), plotx, '+', times(sensor_inds), x_calc, '*');
title('sim-x (+) and x-calc (*) vs times');

subplot(3, 2, 3)
plot(times(sensor_inds), plotz, '+', times(sensor_inds), z_calc, '*');
title('sim-z (+) and z-calc (*) vs times');

subplot(3, 2, 5)
plot(times(sensor_inds), plottheta, '+', times(sensor_inds), theta_calc, '*');
title('sim-theta (+) and theta-calc (*) vs times');

subplot(3, 2, 2)
plot(times(sensor_inds), plotxdot, '+', times(sensor_inds), xdot_calc, '*');
title('sim-xdot (+) and xdot-calc (*) vs times');

subplot(3, 2, 4)
plot(times(sensor_inds), plotzdot, '+', times(sensor_inds), zdot_calc, '*');
title('sim-zdot (+) and zdot-calc (*) vs times');

subplot(3, 2, 6)
plot(times(sensor_inds), plotthetadot, '+', times(sensor_inds), thetadot_calc, '*');
title('sim-thetadot (+) and thetadot-calc (*) vs times');

drawnow;

display(['Cost: ', num2str(prog.objective(zz))]);