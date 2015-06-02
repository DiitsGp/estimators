function [r,xtraj,utraj,ltraj,z,F,info,prog] = contactBasedStateEstimator(xtraj,utraj,ltraj,ljltraj)
%ContactBasedStateEstimator
%   Estimate states given linear acceleration data and angular velocity
%   data otbained from an onboard IMU in a planar falling brick

%% load data

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
% resize them here
x_filt = x_filt(89:end);
y_filt = y_filt(89:end);
theta_filt = theta_filt(89:end);

xdot = xdot(88:end);
ydot = ydot(88:end);
thetadot = thetadot(88:end);

xddot = xddot(87:end);
yddot = yddot(87:end);

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
%    theta_filt(i) = theta_filt(i) + Gdir;
end

%% setup
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.01;
w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');

N = 15;
tf = times(end) - times(1);
inds = round(linspace(1, size(times, 1), N));

x_sensor = [xddot(inds), yddot(inds), thetadot(inds)];
x_gt = [x_filt(inds), y_filt(inds), theta_filt(inds)];

x_calc = zeros(N, 1);
y_calc = zeros(N, 1);
theta_calc = zeros(N, 1);

x_err = zeros(N, 1);
y_err = zeros(N, 1);
theta_err = zeros(N, 1);

urdf = fullfile('CBSE_Window.URDF');
urdf = 'FallingBrick.urdf';
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt);

options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

%% playback original data
v = r.constructVisualizer;
v.display_dt = 0.01;
poses = zeros(6, length(inds));
poses(1:3, :) = [x_filt(inds), y_filt(inds), theta_filt(inds)]';
poses(2, :) = poses(2, :) + 0.05;
dttimes = linspace(times(1), times(end), length(times(inds)));
xtraj_constructed = DTTrajectory(dttimes, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
v.playback(xtraj_constructed);

%% noisy seed
urdf = fullfile('CBSE_Window.URDF');
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt);

x0 = [x_filt(inds(1)); y_filt(inds(1)); theta_filt(inds(1)); xdot(1); ydot(1); thetadot(1)];

x0min = Point(r.getStateFrame());
x0max = Point(r.getStateFrame());

bd = 0.2;
x0min.base_x = x0(1);
x0min.base_z = max(0, x0(2)-bd);
x0min.base_relative_pitch = x0(3) - bd;
x0min.base_xdot = x0(4)-bd;
x0min.base_zdot = x0(5) - bd;
x0min.base_relative_pitchdot = x0(6) - bd;
x0max.base_x = x0(1);
x0max.base_z = x0(2) + bd;
x0max.base_relative_pitch = x0(3) + bd;
x0max.base_xdot = x0(4) + bd;
x0max.base_zdot = x0(5) + bd;
x0max.base_relative_pitchdot = x0(6) + bd;

%% do trajectory optimization
%scale_sequence = [1;.001;0];
prog = ContactImplicitTrajectoryOptimization(r.getManipulator,2,tf,options);
prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
prog = prog.setSolverOptions('snopt','IterationsLimit',200000);

%prog = addStateConstraint(prog, ConstantConstraint(x0),1);
prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);

traj_init.x = PPTrajectory(foh([0, tf/N],[x0,x0]));
[xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);

for i=2:N %length(scale_sequence)
    %  scale = scale_sequence(i);
    display(i);
    options.compl_slack = 0.01;%scale*.01;
    options.lincompl_slack = 0.001;%scale*.001;
    options.jlcompl_slack = 0.01;%scale*.01;
    
    prog = ContactImplicitTrajectoryOptimization(r.getManipulator,i,tf,options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
    prog = prog.setSolverOptions('snopt', 'MajorFeasibilityTolerance', 1e-2);
    prog = prog.setSolverOptions('snopt', 'MajorOptimalityTolerance', 1e-2);
    prog = prog.setSolverOptions('snopt', 'MinorFeasibilityTolerance', 1e-2);
    
    for j = 2:i
        uIMU = [0; 0; 0];
%         uTheta = [0];
        for k = inds(j-1)+1:inds(j)
            dt = times(k) - times(k-1);
            
            uIMU(1) = uIMU(1) + dt*round(xddot(k), 2);
            uIMU(2) = uIMU(2) + dt*round(yddot(k), 2);

%             uTheta = uTheta + dt*round(thetadot(k), 2);
        end
        
        uIMU(3) = round(thetadot(inds(j)), 2);

        IMU_fun = @(x, oldx) IMUcost(x, oldx, uIMU);
        IMUerr_cost = FunctionHandleObjective(2,IMU_fun);
        prog = addCost(prog,IMUerr_cost,{prog.x_inds(4:6, j); prog.x_inds(4:6, j-1)});
        
%         Theta_fun = @(x, oldx) Thetacost(x, oldx, uTheta);
%         Theta_err_cost = FunctionHandleObjective(1, Theta_fun);
%         prog = addCost(prog, Theta_err_cost, {prog.x_inds(3, j); prog.x_inds(3, j-1)});
    end
    
    % initial conditions constraint
    traj_init.x = xtraj;
    traj_init.l = ltraj;
    
%    prog = addStateConstraint(prog, ConstantConstraint(x0),1);
%    prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);

    tic
    [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
    toc
end

%% playback calculated trajectory
urdf = 'FallingBrick.urdf';
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt);
v = r.constructVisualizer;
v.display_dt = 0.01;
v = r.constructVisualizer;
v.display_dt = 0.01;
traj = xtraj.eval(xtraj.getBreaks());
poses = zeros(6, length(inds));
traj(2, :) = traj(2, :)+0.05;
poses(1:3, :) = [traj(1, :)', traj(2, :)', traj(3, :)']';
xtraj_constructed = DTTrajectory(dttimes, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
v.playback(xtraj_constructed, struct('slider', true));

%% plot the final errors
for i = 1:N
    t = times(inds(i));
    xtrajfinal = xtraj.eval(t);
    
    x_calc(i) = xtrajfinal(1);
    y_calc(i) = xtrajfinal(2);
    theta_calc(i) = xtrajfinal(3);
    
    x_err(i) = (x_gt(i, 1) - x_calc(i))/x_gt(i, 1);
    y_err(i) = (x_gt(i, 2) - y_calc(i))/x_gt(i, 2);
    theta_err(i) = (x_gt(i, 3) - theta_calc(i))/x_gt(i, 3);
end

figure
plot(times(inds), x_gt(:, 1), '+', times(inds), x_calc, '*');
title('gt-x (+) and x-calc (*) vs times');

figure
plot(times(inds), x_gt(:, 2), '+', times(inds), y_calc, '*');
title('gt-y (+) and y-calc (*) vs times');

figure
plot(times(inds), x_gt(:, 3), '+', times(inds), theta_calc, '*');
title('gt-theta (+) and theta-calc (*) vs times');

% figure
% stem(times(inds), x_calc, '*');
% title('x-calc at knot points');
% 
% figure
% stem(times(inds), y_calc);
% title('y-calc at knot points');
% 
% figure
% stem(times(inds), theta_calc);
% title('theta-calc at knot points');

drawnow;
%keyboard;

%% cost fun!
    function [f, df] = IMUcost(x, oldx, u)
        costmat = [x(1)-oldx(1); x(2)-oldx(2); x(3)];
        Q = [10/pi, 0, 0; 0, 10/pi, 0; 0, 0, 1]; % Q is a scaling matrix
        f = (costmat-u)'*Q*(costmat-u);
        gradcostmat = [1, 0, 0, -1, 0, 0;...
            0, 1, 0, 0, -1, 0;...
            0, 0, 1, 0, 0, 0];
        df = 2*(costmat-u)'*Q*gradcostmat;
    end

%     function [f, df] = Thetacost(x, oldx, u)
%        costmat = [x-oldx];
%        f = (costmat-u)^2;
%        df = 2*(costmat-u)*[1, -1];
%     end

end

