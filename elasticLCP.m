function [states, times] = elasticLCP(r, x0, tf)
% Solves an elastic LCP as described in Anitescu and Potra 1997

% Re-written to use the forumation in Anitescu 2003, "A Fixed Time-Step
% Approach for Multibody Dynamics with Contact and Friction"

states = x0;
times = 0;
phi_tol = 1e-6;
h = 0.001;
t = 0;
q = x0(1:3);
v = x0(4:6);
vlast = 0;
eps = 0; % coefficient of restitution

while(t < tf)
    [M, C] = r.manipulatorDynamics(q, v, 0);
    kinsol = r.doKinematics(q, v);
    [phi, ~, ~, ~, ~, ~, ~, mu, n, D] = r.contactConstraints(kinsol, false);
    
    contact_inds = find(phi < phi_tol);
    D = vertcat(D{:});
    % Anitescu says it should be Delta = 1/h*phi(contact_inds); but for
    % some reason this induces an energy gain
    Delta = zeros(numel(contact_inds), 1); 
    vn = n(contact_inds, :)*vlast;
    if (isempty(vn))
        Lambda = [];
    else
        Lambda = eps*vn.*ones(numel(contact_inds), 1);
    end
    
    impulse = -h*C; % h*(B*u-C) --> presently no input
    z = anitescuLCP(M, v, impulse, Delta, Lambda, mu(contact_inds),...
        n(contact_inds, :), D([contact_inds, contact_inds+4], :), numel(contact_inds));

    vlast = v;
    v = z(1:3);
    q = q + h*v;
    t = t+h;
    
    states = [states, [q; v]];
    times = [times, t];
end

xtraj_elastic = PPTrajectory(foh(times, states));
vv = r.constructVisualizer();
xtraj_elastic = xtraj_elastic.setOutputFrame(r.getStateFrame);
vv.playback(xtraj_elastic, struct('slider', true));

keyboard;
end


function [z] = anitescuLCP(M, v, impulse, Delta, Lambda, mu, n, D, nC)
if (nC > 0) % to avoid path complaining about empty matrices in the falling case
    display('Collision!');
    mcpmat = [M, -n', -D', zeros(3, nC);...
        n, zeros(nC, 4*nC);
        D, zeros(2*nC, 3*nC), [eye(nC); eye(nC)];...
        zeros(nC, 3), diag(mu), -eye(nC), -eye(nC), zeros(nC, nC)];
    mcpvec = [-M*v-impulse; Delta+Lambda; zeros(3*nC, 1)];
    
    M = mcpmat(1:3, 1:3);
    H = -mcpmat(1:3, 4:end);
    N = mcpmat(4:end, 4:end);
    
    b = mcpvec(1:3);
    w = mcpvec(4:end);
    
    lcpmat = H'*(M\H) + N;
    lcpvec = w-H'*(M\b);
    
    z = pathlcp(lcpmat, lcpvec, zeros(numel(lcpvec),1), inf(numel(lcpvec),1));
    
    x = M\(H*z - b);
    z = [x; z];
else
    lcpmat = M;
    lcpvec = -M*v-impulse;
    
    z = pathlcp(lcpmat, lcpvec, -inf(numel(lcpvec), 1), inf(numel(lcpvec), 1));
end

end