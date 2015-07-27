function [states, times] = elasticLCP(r, x0, tf)
% Solves an elastic LCP as described in Anitescu and Potra 1997

states = x0;
times = 0;
phi_tol = 1e-4;
h = 0.001;
t = 0;
q = x0(1:3);
v = x0(4:6);
eps = 0.1;%0.005; % coefficient of restitution
while(t < tf)
%     if (t > 0.5092)
%         keyboard;
%     end
    [M, C] = r.manipulatorDynamics(q, v, 0);
    
    qnew = q + h*v;
    kinsol = r.doKinematics(qnew);
    [phi, ~, ~, ~, ~, ~, ~, mu, n, D] = r.contactConstraints(kinsol, false);
    D = vertcat(D{:});

    xd = r.dynamics(0, [q;v], 0);
    vnew = v + h*xd(4:6);
    %z = anitescuLCP(M, v, -h*C, mu, n, D, ); % to find vnew
    %vnew = z(1:3)
    
    if ((~ any(phi < phi_tol)) | (v == 0))
        q = qnew;
        v = vnew;
%         disp('no');
        t = t+h;
        states = [states, [q; v]];
        times = [times, t];
    else
%         if (t >  0.440380730094043)
%             keyboard;
%         end
        % linesearch to find collision time with a gap function tolerance
        % of phi_tol
        hhl = 0;
        hhu = h;
        hhm = h/2;
        err = min(phi);

        while (err < 0 || err > phi_tol)
            hhm = (hhl+hhu)/2;
            % added the following block to aid convergence...need to find
            % the sweet spot for the step limit of hhm
%             if (hhm < 1e-5)
%                 hhm = 1e-5;
%                 break;
%             end
            % back to normal operation here
            kinsol = r.doKinematics(q + hhm*v);
            phi = r.contactConstraints(kinsol, true);
            if (min(phi) >= 0)
                hhl = hhm;
            else
                hhu = hhm;
            end
            
            err = min(phi);
        end
        
        hh = hhm; %% hh is the (smaller) timestep to add to t
        qnew = q + hh*v;
        xd = r.dynamics(0, [q; v]); % +++TK: YOU WERE CALLING r.dynamics(0, q, v) here, while dynamics expects t, x, u as arguments
        vminus = v + hh*xd(4:6); % euler step to approximately guess vminus with new timestep
        
%         keyboard;
        
        contact_inds = find(phi < phi_tol);
        z = anitescuLCP(M, vminus, 0, mu(contact_inds), n(contact_inds, :),...
            D([contact_inds, contact_inds+4], :), numel(contact_inds));% to find vc
        vc = z(1:3);
        cn = z(4:3+numel(contact_inds)); % contact normal force
        Fr = eps*n(contact_inds, :)'*cn; % restitution impulse
        
        %%% if it doesn't work, update mu, n, D with contactConstraints at
        %%% qcollision
%         kinsol = r.doKinematics(qnew);
%         [~, ~, ~, ~, ~, ~, ~, mu, n, D] = r.contactConstraints(kinsol, false);
%         D = vertcat(D{:});
        %%% ok
        z = anitescuLCP(M, vc, Fr, mu(contact_inds), n(contact_inds, :),...
            D([contact_inds, contact_inds+4], :), numel(contact_inds));
        vplus = z(1:3);
        v = vplus;
        q = qnew;
        t = t+hh
        states = [states, [q; v]];
        times = [times, t];
%         keyboard;
    end
end

% keyboard;
end


function [z] = anitescuLCP(M, v, impulse, mu, n, D, nC)
    mcpmat = [M, -n', -D', zeros(3, nC);...
        n, zeros(nC, 4*nC);
        D, zeros(2*nC, 3*nC), [eye(nC); eye(nC)];...
        zeros(nC, 3), diag(mu), -eye(nC), -eye(nC), zeros(nC, nC)];
    mcpvec = [-M*v-impulse; zeros(4, 1)];

    M = mcpmat(1:3, 1:3);
    H = -mcpmat(1:3, 4:end);
    N = mcpmat(4:end, 4:end);
    
    b = mcpvec(1:3);

    lcpmat = H'*(M\H) + N;
    lcpvec = -H'*(M\b);
    z = pathlcp(lcpmat, lcpvec, zeros(numel(lcpvec),1), inf(numel(lcpvec),1));
    
    x = M\(H*z - b);
    z = [x; z];
end