options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.0083;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('CBSE_Window.URDF');

G = -5.5;
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt);

tf = 3;
x0 = [0; 1; 5*3.31415926/6; 0; 0; 0];

xtraj_ts = simulate(r, [0 tf], x0);

traj = xtraj_ts.eval(xtraj_ts.getBreaks());
x = traj(1, :);
y = traj(2, :);
theta = traj(3, :);
xdot = traj(4, :);
ydot = traj(5, :);
thetadot = traj(6, :);

times = linspace(0, tf, length(x));
dt = diff(times);
xddot = diff(xdot)./dt(1:end);
yddot = diff(ydot)./dt(1:end);

times = times(1:numel(times)-1) - times(1);
times = times';

xddot = awgn(xddot, 45);
yddot = awgn(yddot, 45);
noisy_thetadot = awgn(thetadot, 45);

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

figure
plot(times(inds), theta(inds), '*');

[r, xtraj, info] = contactBasedStateEstimator(times, x0, xddot, yddot, noisy_thetadot, G, 1);

for i = 1:15
    t = times(inds(i));
    xtrajloop = xtraj.eval(t);
    
    theta_calc(i) = xtrajloop(3);
    xdot_calc(i) = xtrajloop(4);
    ydot_calc(i) = xtrajloop(5);
end

% only plot the results of integration here
figure
plot(times(inds), xdot(inds), '+', times(inds), xdot_calc, '*');
title('sim-xdot (+) and xdot-calc (*) vs times :: scale = 1');

figure
plot(times(inds), ydot(inds), '+', times(inds), ydot_calc, '*');
title('sim-ydot (+) and ydot-calc (*) vs times :: scale = 1');

figure
plot(times(inds), theta(inds), '+', times(inds), theta_calc, '*');
title('sim-theta (+) and theta-calc (*) vs times :: scale = 1');

drawnow;