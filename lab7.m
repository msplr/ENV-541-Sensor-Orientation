% generate reference trajectory
w = pi/100; % [rad]
r = 25; % [m]
t_end = 200; % [s]
dt = 2; % [s]
N = round(t_end/dt + 1);

t = (0:N-1)'*dt;
a = t*w + pi/2;
v = (w*r*[cos(a) sin(a)]);
p = (r*[sin(a) -cos(a)]);

std_gps = [0.5 0.5]; % [m] GPS standard deviation
z_gps = p + std_gps.*randn(N,2);

fprintf('Empirical std dev GPS:%.4f\n', sqrt(sum(var(z_gps-p))))

%% Kalman Filter
Phi = [1 0 dt 0; ...
       0 1 0 dt; ...
       0 0 1 0; ...
       0 0 0 1];
H = [eye(2), zeros(2)];
qvn = 0.05^2; % [(m^2/s^4/Hz]
qve = 0.05^2; % [(m^2/s^4/Hz]
Q = [1/3*dt^3*qvn  0 1/2*dt^2*qvn 0; ...
     0 1/3*dt^3*qve  0 1/2*dt^2*qve; ...
     1/2*dt^2*qvn  0 dt*qvn 0; ...
     0 1/2*dt^2*qve  0 dt*qve];
R = diag(std_gps.^2);

% vectors tracking evolution
x_KF_pred = zeros(N,4);
x_KF = zeros(N,4);
pos_std_dev = zeros(N,1);
innovation = zeros(N,2);

x = [p(1,:) v(1,:)]'; % x = [px py vx vy]'
P = diag([10^2 10^2 0.1^2 0.1^2]);
for i = 1:N
   [xp, Pp] = KF_predict(x,P,Phi,Q);
   z = z_gps(i,:)';
   [x, P] = KF_correct(xp,Pp,z,H,R);

   x_KF_pred(i,:) = xp';
   x_KF(i,:) = x';
   pos_std_dev(i) = sqrt(Pp(1,1) + Pp(2,2)); % KF-predicted positioning quality
   innovation(i,:) = z - H*xp;
end

fprintf('Empirical std dev KF filter:%.4f\n', sqrt(sum(var(xc(:,1:2)-p))))

%% Plots
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

figure
plot(t, pos_std_dev)
title('KF-predicted positioning quality')
xlabel('time [s]'); ylabel('\sigma^{KFp}_{xy} [m]')

figure
plot(t, innovation(:,1)); hold on
plot(t, innovation(:,2))
title('Innovation sequence')
legend('innovation p_{north}', 'innovation p_{east}')
xlabel('time [s]'); ylabel('\sigma^{KFp}_{xy} [m]')

% figure;
% plot(p(:,2), p(:,1))
% hold on
% plot(z_gps(:,2), z_gps(:,1))
% axis equal

% figure;
% plot(p(:,2), p(:,1)); hold on;
% plot(x_KF(:,2), x_KF(:,1))
% plot(x_KF_pred(:,2), x_KF_pred(:,1))
% plot(z_gps(:,2), z_gps(:,1), 'x')
% title('Trajectory'); xlabel('x2 [m]'); ylabel('x1 [m]')
% legend('true position', 'corrected position', 'predicted position', 'GPS measurement')
% axis equal

%% functions
function [x, P] = KF_predict(x,P,Phi,Q)
    x = Phi*x;
    P = Phi*P*Phi' + Q;
end

function [x, P] = KF_correct(x,P,z,H,R)
    K = P*H'*(H*P*H' + R)^-1;
    x = x + K*(z - H*x);
    P = (eye(4) - K*H)*P;
end


