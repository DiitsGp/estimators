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
SNR = getSNR() - 10;
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

%% setup
options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.01;
options.selfCollisions = false;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');

N = 15;
tf = times(end) - times(1);
inds = round(linspace(1, size(times, 1), N));

x_sensor = [xdot(inds), ydot(inds), thetadot(inds)];
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
scale_sequence = [1; 0.5; 0.1; 0];

for scseq=1:length(scale_sequence)
    num = scale_sequence(scseq);
    display(num);
    
    prog = ContactImplicitTrajectoryOptimization(r.getManipulator,2,tf,options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
    
    %prog = addStateConstraint(prog, ConstantConstraint(x0),1);
    prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);
    
    traj_init.x = PPTrajectory(foh([0, tf/N],[x0,x0]));
    [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
    
    for i=2:N %length(scale_sequence)
        scale = scale_sequence(scseq);
        display(i);
        options.compl_slack = scale*.01;
        options.lincompl_slack = scale*.001;
        options.jlcompl_slack = scale*.01;
        
        prog = ContactImplicitTrajectoryOptimization(r.getManipulator,i,tf,options);
        prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
        prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
        prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
        prog = prog.setSolverOptions('snopt', 'FunctionPrecision', 1e-12);
        prog = prog.setSolverOptions('snopt', 'MajorOptimalityTolerance', 5e-4);
        prog = prog.setSolverOptions('snopt', 'MajorFeasibilityTolerance', 1e-4);
        prog = prog.setSolverOptions('snopt', 'MinorFeasibilityTolerance', 1e-4);
        
        for j = 2:i
            uIMU = [0; 0; 0];
            uTHETA = [0];
            for k = inds(j-1)+1:inds(j)
                dt = times(k) - times(k-1);
                
                uIMU(1) = uIMU(1) + dt*round(xddot(k), 2);
                uIMU(2) = uIMU(2) + dt*round(yddot(k), 2);
                
                uTHETA = uTHETA + dt*round(thetadot(k), 2);
            end
            
            uIMU(3) = round(thetadot(inds(j)), 2);
            
            IMU_fun = @(x, oldx) IMUcost(x, oldx, uIMU);
            IMUerr_cost = FunctionHandleObjective(2,IMU_fun);
            prog = addCost(prog,IMUerr_cost,{prog.x_inds(4:6, j); prog.x_inds(4:6, j-1)});
            
            Theta_fun = @(x, oldx) THETAcost(x, oldx, uTHETA);
            Theta_err_cost = FunctionHandleObjective(1, Theta_fun);
            prog = addCost(prog, Theta_err_cost, {prog.x_inds(3, j); prog.x_inds(3, j-1)});
        end
        
        % initial conditions constraint
        traj_init.x = xtraj;
        traj_init.l = ltraj;
        
        %    prog = addStateConstraint(prog, ConstantConstraint(x0),1);
        prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);
        
        tic
        [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
        toc
    end
    
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
    
    figure
    plot(times(inds), x_sensor(:, 1), '+', times(inds), xdot_calc, '*');
    title(['gt-xdot (+) and xdot-calc (*) vs times :: scale = ',num2str(num)]);
    
    figure
    plot(times(inds), x_sensor(:, 2), '+', times(inds), ydot_calc, '*');
    title(['gt-ydot (+) and ydot-calc (*) vs times :: scale = ',num2str(num)]);
    
    figure
    plot(times(inds), x_gt(:, 3), '+', times(inds), theta_calc, '*');
    title(['gt-theta (+) and theta-calc (*) vs times :: scale = ',num2str(num)]);
    
    drawnow;
    
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

    function [f, df] = THETAcost(x, oldx, u)
        costmat = [x-oldx];
        f = 2*(costmat-u)^2;
        df = 2*2*(costmat-u)*[1, -1];
    end

end