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
    theta_filt(i) = theta_filt(i) + Gdir;
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

x0 = [x_filt(inds(1)); y_filt(inds(1)); theta_filt(inds(1)); 0; 0; 0];

urdf = fullfile('CBSE_Window.URDF');
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt);

options = struct();
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

%% do trajectory optimization
%scale_sequence = [1;.001;0];
prog = ContactImplicitTrajectoryOptimization(r.getManipulator,2,tf,options);
prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
prog = prog.setSolverOptions('snopt','IterationsLimit',200000);


prog = addStateConstraint(prog, ConstantConstraint(x0),1);

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
        for k = inds(j-1)+1:inds(j)
            dt = times(k) - times(k-1);
            
            uIMU(1) = uIMU(1) + dt*round(xddot(k), 2);
            uIMU(2) = uIMU(2) + dt*round(yddot(k), 2);
            
        end
        
        uIMU(3) = round(thetadot(inds(j)), 2);
        
        IMU_fun = @(x, oldx) IMUcost(x, oldx, uIMU);
        IMUerr_cost = FunctionHandleObjective(2,IMU_fun);
        prog = addCost(prog,IMUerr_cost,{prog.x_inds(4:6, j); prog.x_inds(4:6, j-1)});
        
        %        uGT = [x_gt(j, 1); x_gt(j, 2); x_gt(j, 3)];
        %        GT_fun = @(x) GTcost(x, uGT);
        %        GTerr_cost = FunctionHandleObjective(1, GT_fun);
        %        prog = prog.addCost(GTerr_cost, {prog.x_inds(1:3, j)'});
    end
    
    % initial conditions constraint
    traj_init.x = xtraj;
    traj_init.l = ltraj;
    
%     prog = addStateConstraint(prog, ConstantConstraint(x0),1);

    tic
    [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
    toc
end

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
stem(times(inds), x_gt(:, 1));
title('gt-x vs times');

figure
stem(times(inds), x_gt(:, 2));
title('gt-y vs times');

figure
stem(times(inds), x_gt(:, 3));
title('gt-theta vs times');

figure
stem(times(inds), x_calc);
title('x-calc at knot points');

figure
stem(times(inds), y_calc);
title('y-calc at knot points');

figure
stem(times(inds), theta_calc);
title('theta-calc at knot points');

drawnow;
keyboard

%% cost fun!
    function [f, df] = IMUcost(x, oldx, u)
        diffvel = [x(1)-oldx(1); x(2)-oldx(2); x(3)];
        Q = [2, 0, 0; 0, 2, 0; 0, 0, 1]; % Q is a scaling matrix
        f = (diffvel-u)'*Q*(diffvel-u);
        graddiffvel = [1, 0, 0, -1, 0, 0;...
            0, 1, 0, 0, -1, 0;...
            0, 0, 1, 0, 0, 0];
        df = 2*(diffvel-u)'*Q*graddiffvel;
    end

%     function [f, df] = GTcost(x, u)
%         f = (x-u)'*(x-u);
%         df = 2*(x-u)';
%     end

end

