function [c] = check_cost()

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.dt = 0.001;
options.use_bullet = false;
options.enable_fastqp = false;
options.ignore_self_collisions = true;

w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
urdf = fullfile('Planar_CBSE_Window.URDF');

G = -10;
p = PlanarRigidBodyManipulator(urdf, options);
p = p.setGravity([0; 0; G]);
r = TimeSteppingRigidBodyManipulator(p, options.dt, options);

tf = 1;
x0 = [0; 0.12; pi/3; 0; 0; 0];

xtraj_ts = simulate(r, [0 tf], x0);

traj = xtraj_ts.eval(xtraj_ts.tt);
x = traj(1, :);
z = traj(2, :);
theta = traj(3, :);
xdot = traj(4, :);
zdot = traj(5, :);
thetadot = traj(6, :);

times = linspace(0, tf, length(x));
dt = diff(times);
sensor_inds = round(linspace(1, numel(xtraj_ts.tt), 120*tf));
sensor_inds = sensor_inds(1:50);
ts = times(sensor_inds);
tempxdot = xdot(sensor_inds);
tempzdot = zdot(sensor_inds);
temptimes = times(sensor_inds);
xddot = diff(tempxdot)./diff(temptimes);
zddot = diff(tempzdot)./diff(temptimes);
theta_dot = thetadot(sensor_inds);
sensor_inds = sensor_inds(1:numel(sensor_inds)-1);

for i = 1:numel(xddot)
    j = sensor_inds(i);
    ang = theta(j);
    R = [cos(ang), sin(ang); -sin(ang), cos(ang)];
    val = R*[xddot(i); zddot(i)];
    xddot(i) = val(1);
    zddot(i) = val(2);
end

noisy_thetadot = theta_dot(2:end);

c = 0;

ts = diff(ts);
measimu = [[0, ts.*xddot]; [0, ts.*zddot]; [0, noisy_thetadot]];
measang = [0, ts.*noisy_thetadot];
plot_inds = round(linspace(1, numel(xtraj_ts.tt), 120*tf));
plot_inds = plot_inds(1:50);
for i=2:50
   [f, ~] = IMUcost2b(traj([4:6, 3], plot_inds(i)), traj([4:6, 3], plot_inds(i-1)), measimu(:, i));
   c = c + f;
   [g, ~] = ANGLEcost2(traj(3, plot_inds(i)), traj(3, plot_inds(i-1)), measang(i));
   c = c + g;
end

    function [f, df] = IMUcost2b(x, oldx, meas)
        %x = [xacc; zacc; thetadot; theta]
        theta = x(4); % positive for R_w^b negative for R_b^w
        R = [cos(theta), sin(theta), 0;...
            -sin(theta), cos(theta), 0;...
            0, 0, 1];
        costmat = [x(1)-oldx(1); x(2)-oldx(2); x(3)];
        Q = 0.01*[1, 0, 0;...
            0, 1, 0;...
            0, 0, 1]; % Q is a scaling matrix
        f = (R*costmat-meas)'*Q*(R*costmat-meas);
        jaccostmat = [cos(theta), sin(theta), 0, -costmat(1)*sin(theta)+costmat(2)*cos(theta), -cos(theta), -sin(theta), 0, 0;...
            -sin(theta), cos(theta), 0, -costmat(1)*cos(theta)-costmat(2)*sin(theta), sin(theta), -cos(theta), 0, 0;...
            0, 0, 1, 0, 0, 0, 0, 0];
        df = 2*(R*costmat-meas)'*Q*jaccostmat;
    end

    function [f, df] = ANGLEcost2(x, oldx, meas)
        angdiff = mod(x-oldx, 2*pi);
        if angdiff <= pi
            costmat = [angdiff];
        else
            costmat = [angdiff-2*pi];
        end
        f = 0.01*(costmat-meas)^2;
        df = 0.01*2*(costmat-meas)*[1, -1];
    end

end

