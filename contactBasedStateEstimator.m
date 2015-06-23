function [r,xtraj,utraj,ltraj,z,F,info,prog] = contactBasedStateEstimator(r, N, x0, data, xtraj)
%ContactBasedStateEstimator
%   Estimate states given linear acceleration data and angular velocity
%   data otbained from an onboard IMU in a planar falling brick

%% setup
options.terrain = RigidBodyFlatTerrain();
options.dt = 0.001;
options.floating = true;
options.selfCollisions = false;
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

times = data(:, 1);
xddot = data(:, 2);
yddot = data(:, 3);
zddot = data(:, 4);
rolldot = data(:, 5);
pitchdot = data(:, 6);
yawdot = data(:, 7);

tf = times(end) - times(1);
inds = round(linspace(1, size(times, 1), N));

%% bounds
x0min = Point(r.getStateFrame());
x0max = Point(r.getStateFrame());

if (r.twoD == 1)
    x0min.base_x = x0(1);
    x0min.base_z = 0;
    x0min.base_relative_pitch = -inf;
    x0min.base_xdot = -inf;
    x0min.base_zdot = -inf;
    x0min.base_relative_pitchdot = -inf;
    
    x0max.base_x = x0(1);
    x0max.base_z = inf;
    x0max.base_relative_pitch = inf;
    x0max.base_xdot = inf;
    x0max.base_zdot = inf;
    x0max.base_relative_pitchdot = inf;
else
    x0min.base_x = x0(1);
    x0min.base_y = x0(2);
    x0min.base_z = 0;
    x0min.base_roll = -inf;
    x0min.base_pitch = 0;
    x0min.base_yaw = 0;
    x0min.base_xdot = -inf;
    x0min.base_ydot = 0;
    x0min.base_zdot = -inf;
    x0min.base_rolldot =  -inf;
    x0min.base_pitchdot = 0;
    x0min.base_yawdot = 0;
    
    x0max.base_x = x0(1);
    x0max.base_y = x0(2);
    x0max.base_z = inf;
    x0max.base_roll = inf;
    x0max.base_pitch = 0;
    x0max.base_yaw = 0;
    x0max.base_xdot = inf;
    x0max.base_ydot = 0;
    x0max.base_zdot = inf;
    x0max.base_rolldot = inf;
    x0max.base_pitchdot = 0;
    x0max.base_yawdot = 0;
end

%% do trajectory optimization
% prog = ContactImplicitTrajectoryOptimization(r.getManipulator,N,tf,options);
% prog = prog.setSolverOptions('snopt','MajorIterationsLimit',500);
% prog = prog.setSolverOptions('snopt','MinorIterationsLimit',500000);
% prog = prog.setSolverOptions('snopt','IterationsLimit',500000);
%
% prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);
%
% traj_init.x = xtraj;
% [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
% trajytraj = xtraj.eval(xtraj.getBreaks());
% display(trajytraj);
% scale_sequence = [1, 0.001, 0];
% 
% for i = 1:numel(scale_sequence)
%     scale = scale_sequence(i);
%     display(['Scale: ', num2str(scale)]);
    options.compl_slack = 0;%scale * 0.01;
    options.lincompl_slack = 0;%scale * 0.001;
    
    prog = ContactImplicitTrajectoryOptimization(r.getManipulator,N,tf,options);
    prog = prog.setSolverOptions('snopt','MajorIterationsLimit',750);
    prog = prog.setSolverOptions('snopt','MinorIterationsLimit',750000);
    prog = prog.setSolverOptions('snopt','IterationsLimit',750000);
    %         prog = prog.setSolverOptions('snopt', 'FunctionPrecision', 1e-12);
    prog = prog.setSolverOptions('snopt', 'MajorOptimalityTolerance', 1e-5);
    prog = prog.setSolverOptions('snopt','print','snopt.out');
%     prog = prog.setCheckGrad(true);
    %         prog = prog.setSolverOptions('snopt', 'MajorFeasibilityTolerance', 1e-6);
    %         prog = prog.setSolverOptions('snopt', 'MinorFeasibilityTolerance', 1e-6);
    
    for kp = 2:N
        display(['Added cost function for knot point: ', num2str(kp)]);
        if (r.twoD == 1)
            measIMU = [0; 0; 0];
            measANG = [0];
            for k = inds(kp-1)+1:inds(kp)
                dt = times(k) - times(k-1);
                measIMU(1) = measIMU(1) + dt*xddot(k);
                measIMU(2) = measIMU(2) + dt*zddot(k);
                measANG = measANG + dt*pitchdot(k);
            end
            measIMU(3) = pitchdot(inds(kp));
            
            IMU_fun = @(x, oldx) IMUcost2(x, oldx, measIMU);
            IMUerr_cost = FunctionHandleObjective(6, IMU_fun);
            
            Angular_fun = @(x, oldx) ANGcost2(x, oldx, measANG);
            Angular_err_cost = FunctionHandleObjective(2, Angular_fun);
            
            prog = addCost(prog,IMUerr_cost,{prog.x_inds(4:6, kp); prog.x_inds(4:6, kp-1)});
            prog = addCost(prog, Angular_err_cost, {prog.x_inds(3, kp); prog.x_inds(3, kp-1)});
        else
            measIMU = [0; 0; 0; 0; 0; 0];
            measANG = [0; 0; 0];
            for k = inds(kp-1)+1:inds(kp)
                dt = times(k) - times(k-1);
                measIMU(1) = measIMU(1) + dt*xddot(k);
                measIMU(2) = measIMU(2) + dt*yddot(k);
                measIMU(3) = measIMU(3) + dt*zddot(k);
                measANG(1) = measANG(1) + dt*rolldot(k);
                measANG(2) = measANG(2) + dt*pitchdot(k);
                measANG(3) = measANG(3) + dt*yawdot(k);
            end
            measIMU(4) = rolldot(inds(kp));
            measIMU(5) = pitchdot(inds(kp));
            measIMU(6) = yawdot(inds(kp));
            
            IMU_fun = @(x, oldx) IMUcost3(x, oldx, measIMU);
            IMUerr_cost = FunctionHandleObjective(12, IMU_fun);
            
            Angular_fun = @(x, oldx) ANGcost3(x, oldx, measANG);
            Angular_err_cost = FunctionHandleObjective(6, Angular_fun);
            
            prog = addCost(prog,IMUerr_cost,{prog.x_inds(7:12, kp); prog.x_inds(7:12, kp-1)});
            prog = addCost(prog, Angular_err_cost, {prog.x_inds(4:6, kp); prog.x_inds(4:6, kp-1)});
        end
    end
    
    traj_init.x = xtraj;
%     if i > 1
%         traj_init.l = ltraj;
%     end
    prog = addStateConstraint(prog, BoundingBoxConstraint(double(x0min), double(x0max)), 1);
    
    tic
    [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
    toc
    assert(info == 1);
% end
v = r.constructVisualizer;
v.display_dt = 0.001;
traj = xtraj.eval(xtraj.getBreaks());
poses = zeros(r.getNumStates, length(inds));
if (r.twoD == 1)
    poses(1:3, :) = [traj(1, :)', traj(2, :)', traj(3, :)']';
else
    poses([1, 2, 3, 4, 5, 6], :) = [traj(1, :)', traj(2, :)', traj(3, :)', traj(4, :)', traj(5, :)', traj(6, :)']';
end
dttimes = linspace(times(1), times(end), length(times(inds)));
xtraj_constructed = DTTrajectory(dttimes, poses);
xtraj_constructed = xtraj_constructed.setOutputFrame(v.getInputFrame);
v.playback(xtraj_constructed, struct('slider', true));

%% cost fun!
    function [f, df] = IMUcost2(x, oldx, meas)
        costmat = [x(1)-oldx(1); x(2)-oldx(2); x(3)];
        Q = [10/pi, 0, 0;...
            0, 10/pi, 0;...
            0, 0, 1]; % Q is a scaling matrix
        f = (costmat-meas)'*Q*(costmat-meas);
        jaccostmat = [1, 0, 0, -1, 0, 0;...
            0, 1, 0, 0, -1, 0;...
            0, 0, 1, 0, 0, 0];
        df = 2*(costmat-meas)'*Q*jaccostmat;
    end

    function [f, df] = ANGcost2(x, oldx, meas)
        costmat = [x-oldx];
        f = (costmat-meas)^2;
        df = 2*(costmat-meas)*[1, -1];
    end

    function [f, df] = IMUcost3(x, oldx, meas)
        costmat = [x(1)-oldx(1); x(2)-oldx(2); x(3)-oldx(3); x(4); x(5); x(6)];
        Q = [10/pi, 0, 0, 0, 0, 0;...
            0, 10/pi, 0, 0, 0, 0;...
            0, 0, 10/pi, 0, 0, 0;...
            0, 0, 0, 1, 0, 0;...
            0, 0, 0, 0, 1, 0;...
            0, 0, 0, 0, 0, 1]; % Q is a scaling matrix
        f = (costmat-meas)'*Q*(costmat-meas);
        jaccostmat = [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;...
            0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;...
            0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0;...
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];
        df = 2*(costmat-meas)'*Q*jaccostmat;
    end

    function [f, df] = ANGcost3(x, oldx, meas)
        costmat = [x(1)-oldx(1); x(2)-oldx(2); x(3)-oldx(3)];
        f = (costmat - meas)'*(costmat - meas);
        jaccostmat = [1, 0, 0, -1, 0, 0;...
            0, 1, 0, 0, -1, 0;...
            0, 0, 1, 0, 0, -1];
        df = 2*(costmat-meas)'*jaccostmat;
    end
end