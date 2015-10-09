options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;
options.use_bullet = false;
options.enable_fastqp = false;
options.ignore_self_collisions = true;
% options.use_new_kinsol = true;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('Planar_CBSE_Window.URDF');

G = -10;
p = PlanarRigidBodyManipulator(urdf, options);
%p = RigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt, options);

tf = 2;
x0 = [0; 1; pi/3; 0; 0; 0];

p = p.compile();

options.integration_method = ContactImplicitTrajectoryOptimization.MIXED;
N = 10;

scale_sequence = [1;.001;0];

for i=1:length(scale_sequence)
  scale = scale_sequence(i);

  options.compl_slack = scale*.01;
  options.lincompl_slack = scale*.001;
  options.jlcompl_slack = scale*.01;
  options.restitution = 0;
  
  prog = ContactImplicitTrajectoryOptimzation(r.getManipulator,N,tf,options);
  prog = prog.setSolverOptions('snopt','MajorIterationsLimit',200);
  prog = prog.setSolverOptions('snopt','MinorIterationsLimit',200000);
  prog = prog.setSolverOptions('snopt','IterationsLimit',200000);
  % prog = prog.setCheckGrad(true);
  
%   snprint('snopt.out');
  
  % initial conditions constraint
  prog = addStateConstraint(prog,ConstantConstraint(x0),1);
  
  if i == 1,
    traj_init.x = PPTrajectory(foh([0,tf],[x0,x0]));
  else
    traj_init.x = xtraj;
    traj_init.l = ltraj;
  end
  [xtraj,utraj,ltraj,~,z,F,info] = solveTraj(prog,tf,traj_init);
end

v = r.constructVisualizer();
xtraj = xtraj.setOutputFrame(r.getStateFrame);
v.playback(xtraj, struct('slider', true));