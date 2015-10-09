options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;
options.use_bullet = false;
options.enable_fastqp = false;
options.ignore_self_collisions = true;
options.restitution = 1;
% options.use_new_kinsol = true;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('Planar_CBSE_Window.URDF');

G = -10;
p = PlanarRigidBodyManipulator(urdf, options);
%p = RigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt, options);

tf = 2;
x0 = [-1.8; 1; pi/4; 0; 0; 0];

p = p.compile();

[states, times] = elasticLCP(p, x0, tf, options.restitution);

xtraj_elastic = PPTrajectory(foh(times, states));
v = r.constructVisualizer();
xtraj_elastic = xtraj_elastic.setOutputFrame(r.getStateFrame);
v.playback(xtraj_elastic, struct('slider', true));
return;

options = struct();
options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;

scale_sequence = [1; 0.001; 0];
for i=1:length(scale_sequence)
  scale = scale_sequence(i);

  options.compl_slack = scale*.01;
  options.lincompl_slack = scale*.001;
  options.jlcompl_slack = scale*.01;
  options.restitution = 1;
  
  prog = ContactImplicitTrajectoryOptimization(p,10,tf,options);
  prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
  prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
  prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
%   prog = prog.setSolverOptions('snopt', 'print', 'snopt.out');
%   prog = prog.setCheckGrad(true);

  % initial conditions constraint
  prog = addStateConstraint(prog,ConstantConstraint(x0),1);
  
  if i == 1,
    traj_init.x = xtraj_elastic;
  else
    traj_init.x = xtraj;
    traj_init.l = ltraj;
  end
  tic
  [xtraj,~,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
  toc
end

xtraj = xtraj.setOutputFrame(r.getStateFrame);
v.display_dt = 0.001;
v.playback(xtraj, struct('slider', true));
