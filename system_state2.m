function ret = system_state2(t, x, u, H, C, g, qd, qd_dot, qd_ddot, s)
    x1 = [x(1); x(3)]; 
    x2 = [x(2); x(4)]; 

    s_value = s(x1, x2, qd(t), qd_dot(t));
    
    dx = zeros(2, 2);

    dx(1, :) = x2';
    dx(2, :) = (H(x1(1), x1(2)) \ u(x1, x2, qd(t), qd_dot(t), qd_ddot(t), s_value) - ...
                H(x1(1), x1(2)) \ C(x1(1), x1(2), x2(1), x2(2)) * x2 - ...
                H(x1(1), x1(2)) \ g(x1(1), x1(2)))';
    
    ret = [dx(1, 1); dx(2, 1); dx(1, 2); dx(2, 2)];
end