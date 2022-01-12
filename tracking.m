clear;
format shortG;

%% parameters

% constant parameters
m1 = 6;
m2 = 4;
l1 = 0.5;
l2 = 0.4;
grav = 9.81;

% system simulation parameters
lc1 = 0.2;
lc2 = 0.1;
I1 = 0.43;
I2 = 0.05;
ml = 0.5;

% max estimation error of parameters
lc1_error = 0.3;
lc2_error = 0.25;
I1_error = 0.48;
I2_error = 0.14;
ml_error = 2;

%% system matrices

% H matrix: inertia matrix of the robot
h11 = @(q1, q2) m1 * lc1^2 + m2 * (lc2 ^ 2 + l1^2 + 2 * l1 * l2 * cos(q2)) + ml * (l2^2 + l1^2 +2 * l1 * l2 * cos(q2)) + I1 + I2;
h12 = @(q1, q2) m2 * lc2 * (lc2 + l1 * cos(q2)) + ml * l2 * (l2 + l1 * cos(q2)) + I2;
h22 = @(q1, q2) lc2^2 * m2 + l2^2 * ml + I2;

H = @(q1, q2) [h11(q1, q2) h12(q1, q2); h12(q1, q2) h22(q1, q2)];

% C matrix: centripetal & Coriolis forces
c11 = @(q1, q2, q1_dot, q2_dot) -l1 * (m2 * lc2 + ml * l2) * sin(q2) * q2_dot;
c12 = @(q1, q2, q1_dot, q2_dot) -l1 * (m2 * lc2 + ml * l2) * sin(q2) * (q2_dot + q1_dot);
c21 = @(q1, q2, q1_dot, q2_dot) l1 * (m2 * lc2 + ml * l2) * sin(q2) * q1_dot;
c22 = @(q1, q2, q1_dot, q2_dot) 0;

C = @(q1, q2, q1_dot, q2_dot) [c11(q1, q2, q1_dot, q2_dot) c12(q1, q2, q1_dot, q2_dot); c21(q1, q2, q1_dot, q2_dot) c22(q1, q2, q1_dot, q2_dot)];

% g array: gravitational forces
g1 = @(q1, q2) (m2 * lc2 + ml * l2) * grav * cos(q1 + q2) + (m2 * l1 + ml * l1 + m1 * lc1) * grav * cos(q1);
g2 = @(q1, q2) (m2 * lc2 + ml * l2) * grav * cos(q1 + q2);

g = @(q1, q2) [g1(q1, q2); g2(q1, q2)];

%% controller matrices

% H maximum error
h11_error = @(q1, q2) m1 * lc1_error + m2 * (lc2_error + 2 * l1 * cos(q2) * lc2_error) + ...
            ml_error * (l1^2 + l2^2 + 2 * l1 * l2 * cos(q2)) + I1_error + I2_error;
h12_error = @(q1, q2) m2 * (lc2_error + l1 * cos(q2) * lc2_error) + ...
            l2 * ml_error * (l2 + l1 * cos(q2)) + I2_error;
h22_error = @(q1, q2) ...
              m2 * lc2_error + l2^2 * ml_error + I2_error;

e1_max = @(q1, q2) ...
         [h11_error(q1, q2) h12_error(q1, q2); ...
          h12_error(q1, q2) h22_error(q1, q2)];

% C maximum error
c11 = @(q1, q2, q1_dot, q2_dot) ...
        -q2_dot;
    
c12 = @(q1, q2, q1_dot, q2_dot) ...
        -(q2_dot + q1_dot);
    
c21 = @(q1, q2, q1_dot, q2_dot) ...
        q1_dot;
    
e2_max = @(q1, q2, q1_dot, q2_dot) ...
          (m2 * lc2_error + ml_error * l2) * (l1 * sin(q2)) * ... 
          [c11(q1, q2, q1_dot, q2_dot) c12(q1, q2, q1_dot, q2_dot); ...
           c21(q1, q2, q1_dot, q2_dot) 0];
      
% g maximum error
g1_error = @(q1, q2) ...
            (m2 * lc2_error + ml_error * l2) * grav * cos(q1 + q2) + ... 
            (ml_error * l1 + m1 * lc1_error) * grav * cos(q1);

g2_error = @(q1, q2) ...
            (m2 * lc2_error + ml_error * l2) * grav * cos(q1 + q2);

e3_max = @(q1, q2) [g1_error(q1, q2); 
                    g2_error(q1, q2)];
                
% e, e_dot (e = q - qd)
e = @(q, qd) q - qd;
e_dot = @(q_dot, qd_dot) q_dot - qd_dot;
                
% r
L = [11 0; 0 11];
C_const = 0.5;
r = @(q, q_dot, qd, qd_dot, qd_ddot) ... 
      e1_max(q(1), q(2)) * abs(qd_ddot - L * e_dot(q_dot, qd_dot)) + ... 
      e2_max(q(1), q(2), q_dot(1), q_dot(2)) * abs(q_dot) + e3_max(q(1), q(2)) + C_const;

% parameter estimations
lc1_hat = 0.25;
lc2_hat = 0.15;
I1_hat = 0.2;
I2_hat = 0.1;
ml_hat = 1;

% C apporximation
c11_hat = @(q1, q2, q1_dot, q2_dot) -l1 * (m2 * lc2_hat + ml_hat * l2) * sin(q2) * q2_dot;
c12_hat = @(q1, q2, q1_dot, q2_dot) -l1 * (m2 * lc2_hat + ml_hat * l2) * sin(q2) * (q2_dot + q1_dot);
c21_hat = @(q1, q2, q1_dot, q2_dot) l1 * (m2 * lc2_hat + ml_hat * l2) * sin(q2) * q1_dot;
c22_hat = @(q1, q2, q1_dot, q2_dot) 0;

C_hat = @(q, q_dot) ...
        [c11_hat(q(1), q(2), q_dot(1), q_dot(2)) c12_hat(q(1), q(2), q_dot(1), q_dot(2)); ...
         c21_hat(q(1), q(2), q_dot(1), q_dot(2)) c22_hat(q(1), q(2), q_dot(1), q_dot(2))];
 
% H apporximation
h11_hat = @(q1, q2) m1 * lc1_hat^2 + m2 * (lc2_hat ^ 2 + l1^2 + 2 * l1 * l2 * cos(q2)) + ml_hat * (l2^2 + l1^2 +2 * l1 * l2 * cos(q2)) + I1_hat + I2_hat;
h12_hat = @(q1, q2) m2 * lc2_hat * (lc2_hat + l1 * cos(q2)) + ml_hat * l2 * (l2 + l1 * cos(q2)) + I2_hat;
h22_hat = @(q1, q2) lc2_hat^2 * m2 + l2^2 * ml_hat + I2_hat;

H_hat = @(q) ...
         [h11_hat(q(1), q(2)) h12_hat(q(1), q(2)); ...
          h12_hat(q(1), q(2)) h22_hat(q(1), q(2))];
      
% g apporximation
g1_hat = @(q1, q2) (m2 * lc2_hat + ml_hat * l2) * grav * cos(q1 + q2) + ...
          (m2 * l1 + ml_hat * l1 + m1 * lc1_hat) * grav * cos(q1);
g2_hat = @(q1, q2) (m2 * lc2_hat + ml_hat * l2) * grav * cos(q1 + q2);

g_hat = @(q) [g1_hat(q(1), q(2)); g2_hat(q(1), q(2))];

% sliding surface of the system
s = @(q, q_dot, qd, qd_dot) ...
        e_dot(q_dot, qd_dot) + L * e(q, qd);

%% controller design

x0 = [pi/3 pi/3; 0 0];
t_range = [0, 10];
epsilon = 0.1;

u = @(q, q_dot, qd, qd_dot, qd_ddot, s_cur) ...
        C_hat(q, q_dot) * q_dot + g_hat(q) + H_hat(q) * qd_ddot - ... 
        H_hat(q) * L * e_dot(q_dot, qd_dot) - ... 
        r(q, q_dot, qd, qd_dot, qd_ddot) .* ... 
        [sat(s_cur(1), epsilon); sat(s_cur(2), epsilon)];

%% compute states by solving the differential equations of the system

qd = @(t) [pi/4 + pi/6 * sin(0.2*pi*t); ...
           -pi/3 + pi/3 * cos(0.2*pi*t)];
   
qd_dot = @(t) [pi^2 / 30 * cos(0.2*pi*t); ...
               -pi^2 / 15 * sin(0.2*pi*t)];
    
qd_ddot = @(t) [-pi^3 / 150 * sin(0.2*pi*t); ...
                -pi^3 / 75 * cos(0.2*pi*t)];
   
[t, x] = ode45(@(t, x) system_state2(t, x, u, H, C, g, qd, qd_dot, qd_ddot, s), t_range, x0);

%% plots

q1 = x(:, 1);
q2 = x(:, 3);
q1_dot = x(:, 2);
q2_dot = x(:, 4);

qd1 = @(t) pi/4 + pi/6 * sin(0.2*pi*t);
qd2 = @(t) -pi/3 + pi/3 * cos(0.2*pi*t);

% system state variables plot
figure(1);
plot(t, q1, 'linewidth', 0.8);
hold on;
plot(t, q2, 'linewidth', 0.8);
hold on
fplot(qd1, [0 20], '--', 'linewidth', 0.8)
hold on
fplot(qd2, [0 20], '--', 'linewidth', 0.8)

legend('q1', 'q2', 'qd1', 'qd2');
xlabel('t (sec)') 
ylabel('rad') 

% deviation of q1 from qd, relative to time
e1 = zeros(size(t));
e2 = zeros(size(t));

for i = 1: length(t)
    e1(i) = abs(q1(i) - qd1(t(i)));
    e2(i) = abs(q2(i) - qd2(t(i)));
end

figure(2);
plot(t, e1, 'linewidth', 0.8);

legend('e1');
xlabel('t (sec)') 
ylabel('rad') 

% deviation of q2 from qd, relative to time
figure(3)
plot(t, e2, 'linewidth', 0.8);

legend('e2');
xlabel('t (sec)') 
ylabel('rad') 

% controller values plot
figure(4);
u_values = zeros(length(t), 2);
for i = 1:size(t) 
    s_value = s([q1(i); q2(i)], [q1_dot(i); q2_dot(i)], ...
        [qd1(t(i)); qd2(t(i))], qd_dot(t(i)));
    u_values(i, :) = u([q1(i); q2(i)], [q1_dot(i); q2_dot(i)], ...
        [qd1(t(i)); qd2(t(i))], qd_dot(t(i)), qd_ddot(t(i)), s_value);
end

plot(t, u_values, 'linewidth', 0.8);
legend('u_1', 'u_2');
xlabel('t (sec)') 
ylabel('exerted force') 

% phase plane plot
qd1_vals = zeros(size(t));
qd2_vals = zeros(size(t));

qd1_dot = @(t) pi^2/30 * cos(0.2*pi*t);
qd2_dot = @(t) -pi^2/15 * sin(0.2*pi*t);

qd1_dot_vals = zeros(size(t));
qd2_dot_vals = zeros(size(t));

for i = 1:length(t)
    qd1_vals(i) = qd1(t(i));
    qd2_vals(i) = qd2(t(i));
    qd1_dot_vals(i) = qd1_dot(t(i));
    qd2_dot_vals(i) = qd2_dot(t(i));
end

e1 = e([q1 q2], [qd1_vals qd2_vals]);
e1 = e1(:, 1);
e1_dot = e_dot([q1_dot q2_dot], [qd1_dot_vals qd2_dot_vals]);
e1_dot = e1_dot(:, 1);

e2 = e([q1 q2], [qd1_vals qd2_vals]);
e2 = e2(:, 2);
e2_dot = e_dot([q1_dot q2_dot], [qd1_dot_vals qd2_dot_vals]);
e2_dot = e2_dot(:, 2);

figure(5);
sliding_surface1 = @(x) -L(1,1) * x;
sliding_surface2 = @(x) -L(2,2) * x;
fplot(sliding_surface1, '--', 'linewidth', 0.8);
hold on;
fplot(sliding_surface2, '--', 'linewidth', 0.8);
hold on;
plot(e1, e1_dot, 'linewidth', 0.8); 
hold on; 
plot(e2, e2_dot, 'linewidth', 0.8);

hl = legend('sliding surface', '$\dot{e_1}$', '$\dot{e_2}$');
set(hl, 'Interpreter', 'latex');

