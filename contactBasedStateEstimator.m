function [r,xtraj,utraj,ltraj,z,F,info,prog] = contactBasedStateEstimator(r, N, x0, data, xtraj)
%ContactBasedStateEstimator
%   Estimate states given linear acceleration data and angular velocity
%   data otbained from an onboard IMU in a planar falling brick

%% setup
options.terrain = RigidBodyFlatTerrain();
options.dt = 0.01;
options.floating = true;
options.selfCollisions = false;
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

times = data(:, 1);
xddot = data(:, 2);
yddot = data(:, 3);
thetadot = data(:, 4);

tf = times(end) - times(1);
inds = round(linspace(1, size(times, 1), N));

%% bounds
% urdf = fullfile('CBSE_Window.URDF');
% p = PlanarRigidBodyManipulator(urdf, options);
% p = p.setGravity([0; 0; G]);
% r = TimeSteppingRigidBodyManipulator(p, options.dt);

x0min = Point(r.getStateFrame());
x0max = Point(r.getStateFrame());

x0min.base_x = x0(1);
x0min.base_z = 0;
x0min.base_relative_pitch = -inf;
x0min.base_xdot = -inf;
x0min.base_zdot = -inf;
x0min.base_relative_pitchdot =  -inf;

x0max.base_x = x0(1);
x0max.base_z = inf;
x0max.base_relative_pitch = inf;
x0max.base_xdot = inf;
x0max.base_zdot = inf;
x0max.base_relative_pitchdot = inf;

%% do trajectory optimization
prog = ContactImplicitTrajectoryOptimization(r.getManipulator,N,tf,options);
prog = prog.setSolverOptions('snopt','MajorIterationsLimit',500);
prog = prog.setSolverOptions('snopt','MinorIterationsLimit',500000);
prog = prog.setSolverOptions('snopt','IterationsLimit',500000);

prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);

traj_init.x = xtraj;
[xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);

for i=2:N
    display(i);
    options.compl_slack = 0;
    options.lincompl_slack = 0;
    
    prog = ContactImplicitTrajectoryOptimization(r.getManipulator,N,tf,options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',750);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',750000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',750000);
    %         prog = prog.setSolverOptions('snopt', 'FunctionPrecision', 1e-12);
    %         prog = prog.setSolverOptions('snopt', 'MajorOptimalityTolerance', 1e-6);
    prog = prog.setSolverOptions('snopt','print','snopt.out');
    %         prog = prog.setCheckGrad(true);
    %         prog = prog.setSolverOptions('snopt', 'MajorFeasibilityTolerance', 1e-6);
    %         prog = prog.setSolverOptions('snopt', 'MinorFeasibilityTolerance', 1e-6);
    
    for j = 2:i
        uIMU = [0; 0; 0];
        uTHETA = [0];
        for k = inds(j-1)+1:inds(j)
            dt = times(k) - times(k-1);
            
            uIMU(1) = uIMU(1) + dt*xddot(k);
            uIMU(2) = uIMU(2) + dt*yddot(k);
            
            uTHETA = uTHETA + dt*thetadot(k);
        end
        
        uIMU(3) = thetadot(inds(j));
        
        IMU_fun = @(x, oldx) IMUcost(x, oldx, uIMU);
        IMUerr_cost = FunctionHandleObjective(6,IMU_fun);
        prog = addCost(prog,IMUerr_cost,{prog.x_inds(4:6, j); prog.x_inds(4:6, j-1)});
        
        Theta_fun = @(x, oldx) THETAcost(x, oldx, uTHETA);
        Theta_err_cost = FunctionHandleObjective(2, Theta_fun);
        prog = addCost(prog, Theta_err_cost, {prog.x_inds(3, j); prog.x_inds(3, j-1)});
    end
    
    traj_init.x = xtraj;
    traj_init.l = ltraj;
    
    prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);
    
    tic
    [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
    toc
    if (info ~= 1 && info ~= 3)
        [c, ceq] = prog.nonlinearConstraints(z);
        display(prog.cin_name);
        display('cin_lb');
        display(prog.cin_lb);
        display('cin_ub');
        display(prog.cin_ub);
        display('c');
        display(c);
    end
    assert(info == 1 | info == 3);
end

v = r.constructVisualizer;
v.display_dt = 0.01;
traj = xtraj.eval(xtraj.getBreaks());
poses = zeros(6, length(inds));
traj(2, :) = traj(2, :);
poses(1:3, :) = [traj(1, :)', traj(2, :)', traj(3, :)']';
dttimes = linspace(times(1), times(end), length(times(inds)));
xtraj_constructed = DTTrajectory(dttimes, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
v.playbackAVI(xtraj_constructed, 'xtraj');

%% cost fun!
    function [f, df] = IMUcost(x, oldx, u)
        costmat = [x(1)-oldx(1); x(2)-oldx(2); x(3)];
        Q = [10/pi, 0, 0; 0, 10/pi, 0; 0, 0, 1]; % Q is a scaling matrix
        f = (costmat-u)'*Q*(costmat-u);
        jaccostmat = [1, 0, 0, -1, 0, 0;...
            0, 1, 0, 0, -1, 0;...
            0, 0, 1, 0, 0, 0];
        df = 2*(costmat-u)'*Q*jaccostmat;
    end

    function [f, df] = THETAcost(x, oldx, u)
        costmat = [x-oldx];
        f = 2*(costmat-u)^2;
        df = 2*2*(costmat-u)*[1, -1];
    end

end