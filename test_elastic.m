options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;
options.use_bullet = false;
options.enable_fastqp = false;
options.ignore_self_collisions = true;
options.use_new_kinsol = true;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('Planar_CBSE_Window.URDF');

G = -10;
p = PlanarRigidBodyManipulator(urdf, options);
%p = RigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt, options);

tf = 10;
x0 = [0; 1; pi/4; 0.5; 0; pi];

p = p.compile();

[states, times] = elasticLCP(p, x0, 2);

xtraj_elastic = PPTrajectory(foh(times, states));
v = r.constructVisualizer();
xtraj_elastic = xtraj_elastic.setOutputFrame(r.getStateFrame);
v.playback(xtraj_elastic, struct('slider', true))
