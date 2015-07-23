function [states, times] = elasticLCP(r, x0, tf)
% Solves an elastic LCP as described in Anitescu and Potra 1997

states = x0;
times = 0;
phi_tol = 1e-6;
h = 0.001;
t = 0;
q = x0(1:3);
v = x0(4:6);
eps = 0.5; % coefficient of restitution
while(t < tf)
    [M, C] = r.manipulatorDynamics(q, v, 0);
    
    qnew = q + h*v;
    kinsol = r.doKinematics(qnew);
    [phi, ~, ~, ~, ~, ~, ~, mu, n, D] = r.contactConstraints(kinsol, false);
    D = vertcat(D{:});

    z = anitescuLCP(M, v, -h*C, mu, n, D); % to find vnew
    vnew = z(1:3)
    
    if (~ any(phi <= 0))
        q = qnew;
        v = vnew;
%         disp('no');
        t = t+h
        states = [states, [q; v]];
        times = [times, t];
    else
        % linesearch to find collision time with a gap function tolerance
        % of phi_tol
        hhl = 0;
        hhu = h;
        hhm = h/2;
        err = abs(min(phi));
        while (err > phi_tol)
            hhm = (hhl+hhu)/2;
            % added the following block to aid convergence...need to find
            % the sweet spot for the step limit of hhm
            if (hhm < 1e-8)
                hhm = 1e-8;
                break;
            end
            % back to normal operation here
            kinsol = r.doKinematics(q + hhm*v);
            phi = r.contactConstraints(kinsol, false);
            if (min(phi) >= 0)
                hhl = hhm;
            else
                hhu = hhm;
            end

            err = abs(min(phi));
        end
        
        hh = hhm; %% hh is the (smaller) timestep to add to t
        qnew = q + hh*v;
        xd = r.dynamics(0, [q; v]); % +++TK: YOU WERE CALLING r.dynamics(0, q, v) here, while dynamics expects t, x, u as arguments
        vminus = v + hh*xd(4:6); % euler step to approximately guess vminus with new timestep
        
        z = anitescuLCP(M, vminus, 0, mu, n, D); % to find vc
        vc = z(1:3);
        cn = z(4:7); % contact normal force
        Fr = sum(eps*n'*cn); % restitution impulse
        
        %%% if it doesn't work, update mu, n, D with contactConstraints at
        %%% qcollision
        kinsol = r.doKinematics(qnew);
        [~, ~, ~, ~, ~, ~, ~, mu, n, D] = r.contactConstraints(kinsol, false);
        D = vertcat(D{:});
        %%% ok
        z = anitescuLCP(M, vc, Fr, mu, n, D);
        vplus = z(1:3);
        v = vplus;
        q = qnew;
        t = t+hh
        states = [states, [q; v]];
        times = [times, t];
%         keyboard;
    end
end


end

function [z] = anitescuLCP(M, v, impulse, mu, n, D)
    mcpmat = [M, -n', -D', zeros(3, 4);...
        n, zeros(4, 16);...
        D, zeros(8, 12), [eye(4);eye(4)];...
        zeros(4, 3), diag(mu), -eye(4), -eye(4), zeros(4, 4)];
    mcpvec = [-M*v - impulse; zeros(16, 1)];
    A = mcpmat(1:3, 1:3);
    B = mcpmat(1:3, 4:end);
    C = -B';
    M = mcpmat(4:end, 4:end);
    
    b = mcpvec(1:3);

    lcpmat = C*(A\B) + M;
    lcpvec = C*(A\b);
    
    z = pathlcp(lcpmat, lcpvec, zeros(16,1), inf(16,1));
    
    x = A\(B*z - b);
    
    z = [x; z];
end


