%% load data
format long g;
OptiData = csvread('Mar312015Sample01_2d.csv'); %time, x, y, theta
%checkDependency('lcmgl');
%lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'FallingPlate');

%% filter and turn into xddot, yddot, thetadot
[b, a] = butter(1, 0.2);

times = OptiData(:, 1);
x_filt = filtfilt(b, a, OptiData(:, 2));
y_filt = filtfilt(b, a, OptiData(:, 3));
y_filt = y_filt - y_filt(end);
theta_filt = filtfilt(b, a, OptiData(:, 4))*pi/180;
dt = diff(times);
xdot = diff(x_filt) ./ dt(1:end);
ydot = diff(y_filt) ./ dt(1:end);
xddot = diff(xdot) ./ dt(2:end);
yddot = diff(ydot) ./ dt(2:end);
thetadot = diff(theta_filt) ./ dt;

% because the samples are motionless in the beginning it's ok for me to
%   resize them here

x_filt = x_filt(89:end);
y_filt = y_filt(89:end);
theta_filt = theta_filt(89:end);

xdot = xdot(88:end);
ydot = ydot(88:end);
thetadot = thetadot(88:end);

xddot = xddot(87:end);
yddot = yddot(87:end);

% add white noise
% I let the IMU sit on a table and measured the gravitational acceleration.
%   Using a known value for g = -9.80665, I calculated the SNR of the IMU as
%   50.7912.
SNR = getSNR() - 10; % make it more noisy than the actual sensor just for funsies
thetadot = awgn(thetadot, SNR);
xddot = awgn(xddot, SNR);
yddot = awgn(yddot, SNR);

times = times(89:end);
times = times - times(1);

G = zeros(53, 1);
Gdir = zeros(53, 1);
for i=25:77
    G(i-24) = norm([xddot(i), yddot(i)]);
    Gdir(i-24) = atan2(yddot(i), xddot(i));
end
G = -mean(G);
Gdir = mean(Gdir) + pi;
Rot = [cos(Gdir), -sin(Gdir), 0; sin(Gdir), cos(Gdir), 0; 0, 0, 0];
for i = 1:length(x_filt)
    State = Rot*[x_filt(i); y_filt(i); 0];
    x_filt(i) = State(1);
    y_filt(i) = State(2);
    Vel = Rot*[xdot(i); ydot(i); 0];
    xdot(i) = Vel(1);
    ydot(i) = Vel(2);
    Acc = Rot*[xddot(i); yddot(i); 0];
    xddot(i) = Acc(1);
    yddot(i) = Acc(2);
end

N = 15;
tf = times(end) - times(1);
inds = round(linspace(1, size(times, 1), N));
x0 = [x_filt(inds(1)); y_filt(inds(1)); theta_filt(inds(1)); xdot(1); ydot(1); thetadot(1)];
[r, xtraj, info] = contactBasedStateEstimator(times, x0, xddot, yddot, thetadot, G, 1);

xdot_gt = [xdot(inds), ydot(inds), thetadot(inds)];
x_gt = [x_filt(inds), y_filt(inds), theta_filt(inds)];

x_calc = zeros(N, 1);
y_calc = zeros(N, 1);
theta_calc = zeros(N, 1);

x_err = zeros(N, 1);
y_err = zeros(N, 1);
theta_err = zeros(N, 1);

%% plot the final errors
for i = 1:N
    t = times(inds(i));
    xtrajloop = xtraj.eval(t);
    
    x_calc(i) = xtrajloop(1);
    y_calc(i) = xtrajloop(2);
    theta_calc(i) = xtrajloop(3);
    xdot_calc(i) = xtrajloop(4);
    ydot_calc(i) = xtrajloop(5);
    
    x_err(i) = (x_gt(i, 1) - x_calc(i))/x_gt(i, 1);
    y_err(i) = (x_gt(i, 2) - y_calc(i))/x_gt(i, 2);
    theta_err(i) = (x_gt(i, 3) - theta_calc(i))/x_gt(i, 3);
end

% only plot the results of integration here
figure
plot(times(inds), xdot_gt(:, 1), '+', times(inds), xdot_calc, '*');
title('gt-xdot (+) and xdot-calc (*) vs times :: scale = 0.1');

figure
plot(times(inds), xdot_gt(:, 2), '+', times(inds), ydot_calc, '*');
title('gt-ydot (+) and ydot-calc (*) vs times :: scale = 0.1');

figure
plot(times(inds), x_gt(:, 3), '+', times(inds), theta_calc, '*');
title('gt-theta (+) and theta-calc (*) vs times :: scale = 0.1');

drawnow;