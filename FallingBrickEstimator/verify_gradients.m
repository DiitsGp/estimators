function verify_gradients()
x1 = [randn(4, 1)];
x2 = [randn(4, 1)];
meas = zeros(3, 1);

eps = 1e-6;

[f, df] = IMUcost(x1, x2, meas);
k = size(df);
dg = zeros(k);
for i = 1:2*numel(x1)
    if i <= numel(x1)
        temp = x1;
        temp(i) = x1(i)+eps;
        [g, l] = IMUcost(temp, x2, meas);
        dg(i) = (g-f)/eps;
    else
        temp = x2;
        temp(i-numel(x1)) = x2(i-numel(x1))+eps;
        [g, l] = IMUcost(x1, temp, meas);
        dg(i) = (g-f)/eps;
    end
end
display(df);
display(dg);

    function [f, df] = IMUcost(x, oldx, meas)
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

end