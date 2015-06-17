options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('CBSE_Window.URDF');

G = -4;
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt);

N = 15;
tf = 2.87;
x0 = [0; 1; 2*pi/3; 0; 0; 0];

xtraj_ts = simulate(r, [0 tf], x0);

traj = xtraj_ts.eval(xtraj_ts.getBreaks());
x = traj(1, :);
y = traj(2, :);
theta = traj(3, :) - pi;
xdot = traj(4, :);
ydot = traj(5, :);
thetadot = traj(6, :);

times = linspace(0, tf, length(x));
dt = diff(times);
sensor_inds = round(linspace(1, numel(xdot)-1, round(120*tf))-1)+1;
tempxdot = xdot(sensor_inds);
tempydot = ydot(sensor_inds);
tempdt = dt(sensor_inds)*numel(times)/numel(sensor_inds);
xddot = diff(tempxdot)./tempdt(1:numel(tempdt)-1);
yddot = diff(tempydot)./tempdt(1:numel(tempdt)-1);
sensor_inds = sensor_inds(1:numel(sensor_inds)-1);

times = times(1:numel(times)-1) - times(1);
times = times';

xddot = awgn(xddot, 45);
yddot = awgn(yddot, 45);
noisy_thetadot = awgn(thetadot(sensor_inds), 45);

inds = round(linspace(1, size(times, 1), 15));

v = r.constructVisualizer;
v.display_dt = 0.01;
v = r.constructVisualizer;
poses = zeros(6, numel(x));
poses(1:3, :) = [x; y; theta];
dttimes = linspace(times(1), times(end), numel(x));
xtraj_constructed = DTTrajectory(dttimes, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
v.playback(xtraj_constructed, struct('slider', true));

% figure
% plot(times(inds), theta(inds), '*');
traj_sim = simulate(r, [0 tf], x0);
ts_sim = traj_sim.getBreaks();
traj_sim = PPTrajectory(foh(ts_sim,traj_sim.eval(ts_sim)));
ts = times(sensor_inds);
data = [ts, xddot', yddot', noisy_thetadot'];
[r, xtraj, info] = contactBasedStateEstimator(r, N, x0, data, traj_sim);

for i = 1:N
    t = times(inds(i));
    xtrajloop = xtraj.eval(t);
    
    theta_calc(i) = xtrajloop(3);
    xdot_calc(i) = xtrajloop(4);
    ydot_calc(i) = xtrajloop(5);
end

% only plot the results of integration here
figure
plot(times(inds), xdot(inds), '+', times(inds), xdot_calc, '*');
title('sim-xdot (+) and xdot-calc (*) vs times');

figure
plot(times(inds), ydot(inds), '+', times(inds), ydot_calc, '*');
title('sim-ydot (+) and ydot-calc (*) vs times');

figure
plot(times(inds), theta(inds), '+', times(inds), theta_calc, '*');
title('sim-theta (+) and theta-calc (*) vs times');


drawnow;