function [states, times] = elasticLCP(r, x0, tf, eps)
% Solves an elastic LCP as described below:
% 
% v_{k+1} = v_{k} + H^-1*(B_{k}*u-C_{k} + J_{k}'*(lambda_c + lambda_d))
% [[lambda_d = e*lambda_c where e is restitution coefficient]]
% v_{c} = v_{k} + H^-1*(B_{k}*u-C_{k}+J_{k}'*lambda_c)
% 0 <= lambda_c \perp phi_{k} + h*J_{k}*v_{c} >= 0
% 
eps = 1;
states = x0;
times = 0;
phi_tol = 1e-3;
h = 1e-2;
t = 0;
q = x0(1:3);
v = x0(4:6);

while(t < tf)
    z = elasticLCP(q, v, zeros(3, 1));
    v = z(1:3);
    q = q+h*v;
    t = t+h
    
    states = [states, [q; v]];
    times = [times, t];
end

xtraj_elastic = PPTrajectory(foh(times, states));
vv = r.constructVisualizer();
xtraj_elastic = xtraj_elastic.setOutputFrame(r.getStateFrame);
vv.playback(xtraj_elastic, struct('slider', true));

function [z] = elasticLCP(q, vk, u)
    % nC is total number of (active and inactive) contact points
    [H, C, B] = r.manipulatorDynamics(q, vk, 0);
    kinsol = r.doKinematics(q, vk);
    [phi, ~, ~, ~, ~, ~, ~, mu, n, D] = r.contactConstraints(kinsol, false);
    D = vertcat(D{:});
    
    nC = r.getNumContactPairs;
    
    MLCPmat = [H, [-n', -D']*h, zeros(3, nC);...
        n*h, zeros(nC, 4*nC);...
        D, zeros(2*nC, 3*nC), [eye(nC); eye(nC)];...
        zeros(nC, 3), [diag(mu), -eye(nC), -eye(nC)]*h, zeros(nC, nC)];
    MLCPvec = [-H*vk+h*C; eps*n*vk; zeros(3*nC, 1)];
    
    % hardcoded indices for now, with the 2d case
    H = MLCPmat(1:3, 1:3);
    G = -MLCPmat(1:3, 4:end);
    F = MLCPmat(4:end, 1:3);
    N = MLCPmat(4:end, 4:end);
    
    k = MLCPvec(1:3);
    l = MLCPvec(4:end);
    
    LCPmat = F*(H\G) + N;
    LCPvec =  F*(H\k) - l;

    z = pathlcp(LCPmat, LCPvec, zeros(numel(LCPvec),1), inf(numel(LCPvec),1));
    
    lambda_N = z(1:nC);
    lambda_T = z(1+nC:3*nC);
    vkplus1 = H\(G*z - k); % ideally approximately zero for now
%     delta_v_rest = h*eps*[n', D']*[lambda_N; lambda_T];
%     vkplus1 = vc + H\delta_v_rest; % for now this should be equal to -vk
%     if (vkplus1 ~= vc)
%         display('collision!');
% %         keyboard;
%     end
%     
    z = [vkplus1; z];
end
end