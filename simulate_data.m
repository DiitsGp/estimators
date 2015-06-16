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
xddot = diff(xdot)./dt(1:end);
yddot = diff(ydot)./dt(1:end);

times = times(1:numel(times)-1) - times(1);
times = times';

xddot = awgn(xddot, 45);
yddot = awgn(yddot, 45);
noisy_thetadot = awgn(thetadot, 45);

sensor_inds = round(linspace(1, numel(xddot),  round(120*tf)));
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
xtraj = PPTrajectory(foh([0, tf/N],[x0, x0]));
[r, xtraj, info] = contactBasedStateEstimator(times(sensor_inds), x0, xddot(sensor_inds), yddot(sensor_inds), noisy_thetadot(sensor_inds), G);

for i = 1:15
    t = times(inds(i));
    xtrajloop = xtraj.eval(t);
    
    theta_calc(i) = xtrajloop(3);
    xdot_calc(i) = xtrajloop(4);
    ydot_calc(i) = xtrajloop(5);
end

% only plot the results of integration here
figure
plot(times(inds), xdot(inds), '+');
title('sim-xdot (+)');

figure
plot(times(inds), ydot(inds), '+');
title('sim-ydot (+)');

figure
plot(times(inds), theta(inds), '+');
title('sim-theta (+)');

drawnow;