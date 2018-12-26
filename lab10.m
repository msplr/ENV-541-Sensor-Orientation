% Lab 10
% author: Michael Spieler

ERROR_FEEDBACK = true;
%ERROR_FEEDBACK = false;
rng(1234);

w0 = pi/100; % [rad]
r = 500; % [m]
t_end = round(2*pi/w0); % [s]
f = 100;
dt = 1/f; % [s] Kalman filter update time
t = (dt:dt:t_end)'; % [s] timestamp of KF updates
dt_GPS = 2; % [s]
t_GPS = (dt_GPS:dt_GPS:t_end)'; % [s] timestamp of GPS measurements

% reference trajectory
a_ref = @(t) t*w0 + pi/2;
v_ref = @(t) (w0*r*[cos(a_ref(t)) sin(a_ref(t))]);
p_ref = @(t) (r*[sin(a_ref(t)) -cos(a_ref(t))]);

% Noise parameters
g = 9.81; % [m/s^2]
gyro_std_bias = -400*pi/(180*3600); % [rad/s]
gyro_std_wn = 0.1*pi*sqrt(1/dt)/(180*60); % [rad/s/sample]
gyro_std_GM = 0.01*pi/180;
gyro_beta_GM = 1/30; % [1/s]
acc_std_wn = 5e-5*g*sqrt(1/dt); % [m/s^2/sample]
acc_std_GM = 2e-4*g; % [m/s^2/sqrt(Hz)]
acc_beta_GM = 1/60; % [1/s]
GPS_std_wn = 1; % [m]

% simulate measurements
gyro_bc = gyro_std_bias * randn();
gyro_wn = gyro_std_wn * randn(length(t),1);
gyro_GM = GaussMarkov_1st_order(gyro_std_GM*randn(length(t),1), dt, gyro_beta_GM);
gyro = w0 + gyro_GM + gyro_wn + gyro_bc;

acc_wn = acc_std_wn * randn(length(t),2);
acc_GM = GaussMarkov_1st_order(acc_std_GM*randn(length(t),2), dt, acc_beta_GM);
accel = [0, r*w0^2] + acc_GM + acc_wn;

GPS_wn = GPS_std_wn * randn(length(t_GPS), 2);
GPS = p_ref(t_GPS) + GPS_wn;

% functions for F and G matrices
Rbm = @(a) [cos(a), -sin(a); sin(a), cos(a)];
F11 = @(a,f) [zeros(1,5);
              Rbm(a)*[-f(2); f(1)], zeros(2,4);
              0, 1, 0, 0, 0;
              0, 0, 1, 0, 0];
F12 = @(a) [1, 1, 0, 0;
            zeros(2), Rbm(a);
            zeros(2,4)];
F21 = zeros(4,5);
F22 = diag([0, -gyro_beta_GM, -acc_beta_GM, -acc_beta_GM]);
get_F = @(a,f) [F11(a,f), F12(a);
                F21, F22];

G11 = @(a) [1, 0, 0; zeros(2,1), Rbm(a); zeros(2,3)];
G12 = zeros(5,3); G21 = zeros(4,3); G22 = [zeros(1,3); eye(3)];
get_G = @(a) [G11(a), G12;
              G21, G22];

% initialize KF variables
x1 = [a_ref(0) v_ref(0) p_ref(0)]';
dx = zeros(9,1);
z_gyro = zeros(length(t),1);
z_accel = zeros(length(t),2);
P = diag([2*pi/180, 5, 5, 10, 10, 0.05*pi/180, 0.01*pi/180, 300e-6*g, 300e-6*g].^2);
e = 0;
W = diag([gyro_std_wn^2, acc_std_wn^2, acc_std_wn^2, 2*gyro_std_GM^2*gyro_beta_GM, ...
          2*acc_std_GM^2*acc_beta_GM, 2*acc_std_GM^2*acc_beta_GM]);
R = diag([GPS_std_wn, GPS_std_wn].^2);
H = zeros(2,9); H(1,4) = 1; H(2,5) = 1;

x1_log = zeros(length(t),5);
x1_KF_log = zeros(length(t_GPS),5);
dz_log = zeros(length(t_GPS),2);
dx_log = zeros(length(t),9);
P_log = zeros(length(t),length(dx),length(dx));
i_GPS = 1;
% simulate KF
for i = 1:length(t)
    z_gyro(i) = gyro(i);
    z_accel(i,:) = accel(i,:);

    if ERROR_FEEDBACK
        e = e + dx(6:end);
        dx(6:end) = 0;
        z_gyro(i) = z_gyro(i) - e(1) - e(2);
        z_accel(i,:) = z_accel(i,:) - e(3:4)';
    end

    % motion model using strapdown INS
    x1 = strapdown_INS(x1, z_gyro, z_accel, i, dt);

    F = get_F(x1(1), z_accel(i,:)'); % uses angle a and f1, f2
    G = get_G(x1(1)); % uses angle a
    [Phi, Q] = get_Phi_and_Q(F, G, W, dt);

    % prediction step
    dx = Phi*dx;
    P = Phi*P*Phi' + Q;
    if t(i) >= t_GPS(i_GPS)
        % correction step
        dz = GPS(i_GPS,:)' - x1(4:end);
        K = P*H'*(H*P*H' + R)^-1;
        dx = dx + K*(dz - H*dx);

        I = eye(length(dx));
        P = (I - K*H)*P;

        x1_KF_log(i_GPS,:) = x1 + dx(1:5);
        %x1 = x1 + dx(1:5);
        %dx(1:5) = 0;

        i_GPS = i_GPS + 1;
    end
    
    dx_log(i,:) = dx;

    x1 = x1 + dx(1:5);
    dx(1:5) = 0;

    P_log(i,:,:) = P;
    x1_log(i,:) = x1';
end

fprintf('Empirical std dev GPS: %.4f m\n', sqrt(sum(var(GPS-p_ref(t_GPS)))))
fprintf('Empirical std dev KF filter: %.4f m\n', sqrt(sum(var(x1_log(:,1:2)-p_ref(t)))))
fprintf('Empirical std dev KF filter: %.4f m\n', sqrt(sum(var(x1_KF_log(:,1:2)-p_ref(t_GPS)))))

%%
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

%plot_pos_and_heading(x1_log, t, P_log);
plot_traj(p_ref(t), x1_log(:,4:5), GPS);
figure; plot(dx_log(:,1:5)); legend('da','dv1','dv2','dp1','dp2')
figure; plot(dx_log(:,6:9)); legend('bc','bg','ba1','ba2')
%% functions
function x = strapdown_INS(x, gyro, accel, i, dt)
    a_ = x(1); v_ = x(2:3); p_ = x(4:5);
    Rbm = @(a) [cos(a) -sin(a); sin(a), cos(a)];
    g = gyro(i);
    f = accel(i,:)';
    if i == 1
        g_ = g;
        f_ = f;
    else
        g_ = gyro(i-1);
        f_ = accel(i-1,:)';
    end
    a = a_ + dt*(g_ + g)/2;
    v = v_ + dt*(Rbm(a_)*f_ + Rbm(a)*f)/2;
    p = p_ + dt*(v_ + v)/2;
    x = [a; v; p];
end

function [Phi, Qk] = get_Phi_and_Q(F, G, W, dt)
    n = length(F);
    A = [-F, G*W*G'; ...
              zeros(n,n), F'];
    B = expm(dt*A);
    Phi = B(n+1:end,n+1:end)';
    Qk = Phi*B(1:n,n+1:end);
end

function x = GaussMarkov_1st_order(w, dt, beta)
    x = zeros(size(w));
    xk = 0;
    for i = 1:length(w)
        xk = exp(-beta*dt).*xk + w(i,:);
       x(i,:) = xk;
    end
end

function [] = plot_traj(p, x_KF, z_gps)
    figure;
    plot(p(:,2), p(:,1)); hold on;
    plot(x_KF(:,2), x_KF(:,1))
    plot(z_gps(:,2), z_gps(:,1), 'x')
    title('Trajectory'); xlabel('x2 [m]'); ylabel('x1 [m]')
    legend('true position', 'estimated position', 'GPS measurement')
    axis equal
end

function [] = plot_signal_with_std(t, y, std)
    yl = y - std;
    yu = y + std;
    fill([t; flipud(t)], [yu; flipud(yl)], [.9 .9 .9], 'linestyle', 'none')
    hold on;
    plot(t,y);
end

function [] = plot_pos_and_heading(x, t, P)
    figure;
    subplot(3,1,1)
    plot_signal_with_std(t, 180/pi*x(:,1), sqrt(P(:,1,1)))
%    plot(t,180/pi*x(:,1)); hold on;
%    plot(t,180/pi*(x(:,1)+sqrt(P(:,1,1))));
%    plot(t,180/pi*(x(:,1)-sqrt(P(:,1,1))));
    subplot(3,1,2)
    plot_signal_with_std(t, x(:,2), sqrt(P(:,2,2)));
    plot_signal_with_std(t, x(:,3), sqrt(P(:,3,3)));
    % plot(t,x(:,2)); hold on;
    % plot(t,x(:,2)+sqrt(P(:,2,2)));
    % plot(t,x(:,2)-sqrt(P(:,2,2)));
    % plot(t,x(:,3))
    % plot(t,x(:,3)+sqrt(P(:,3,3)));
    % plot(t,x(:,3)-sqrt(P(:,3,3)));
    subplot(3,1,3)
    plot(t,x(:,4)); hold on;
    plot(t,x(:,4)+sqrt(P(:,4,4)));
    plot(t,x(:,4)-sqrt(P(:,4,4)));
    plot(t,x(:,5))
    plot(t,x(:,5)+sqrt(P(:,5,5)));
    plot(t,x(:,5)-sqrt(P(:,5,5)));
end
