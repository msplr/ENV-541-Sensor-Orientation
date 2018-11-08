% generate reference trajectory
w = pi/100;
r = 500;
% 100Hz
f_int = 100;
t_end = 200;
dt = 1/f_int;
N = t_end*f_int + 1;

a = (0:20000)'*pi/100*0.01 + pi/2;
v = (w*r*[cos(a) sin(a)]);
p = (r*[sin(a) -cos(a)]);

% sensor error flags
gyro_err_random_bias = 0;
gyro_err_gauss_markov = 0;
gyro_err_white_noise = 0;

acc_err_random_bias = 0;
acc_err_white_noise = 0;

rng(42);

% stochastic error values
gyro_std_bias = 10 /180*pi/3600; % [rad/s]
gyro_std_GM = 0.005 *pi/180; % [rad/s/sqrt(Hz)]
gyro_beta_GM = 1/100; % [1/s]
gyro_std_wn = 0.1 *pi/180 *10/60; % [rad/s/sample]

acc_std_bias = 1e-3 *9.807; % [m/s^2]
acc_std_wn = 9.807*50e-6*10; % [m/s^2/sample]

% generate errors
beta = gyro_beta_GM;
sigma = gyro_std_GM;
wn_GM = sqrt((1-exp(-2*beta*dt))*sigma^2) * randn(N,1);
gyro_GM = GMP(wn_GM, dt, beta) * gyro_err_gauss_markov;

gyro_bc = gyro_std_bias * randn() * gyro_err_random_bias;
gyro_wn = gyro_std_wn * randn(N,1) * gyro_err_white_noise;

w0 = pi/100;
w = w0 + gyro_GM + gyro_wn + gyro_bc;

acc_wn = acc_std_wn * randn(N,2) * acc_err_white_noise;
acc_bc = acc_std_bias * randn(1,2) * acc_err_random_bias;

f = [0, r*w0^2] + acc_bc + acc_wn;

%% Error plots
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

% [atest, vtest, ptest] = trapezoidal_int(a(1), v(1,:)', p(1,:)', f, w, dt);
% t = dt * (0:length(p)-1);
% 
% desc = make_description(gyro_err_random_bias, gyro_err_gauss_markov, gyro_err_white_noise, acc_err_random_bias, acc_err_white_noise);
% 
% plot_trajectory('traj', ptest)
% plot_error(desc, atest-a, vtest-v, ptest-p, t)

%%
desc = make_description(gyro_err_random_bias, gyro_err_gauss_markov, gyro_err_white_noise, acc_err_random_bias, acc_err_white_noise);

sim_and_plot(f, w, a, v, p, dt, desc)
printpdf(gcf, strcat('05/',datestr(now, 'yyyymmdd-HHMMSS'),'.pdf'),1,1.3)

%% functions
function [asim, vsim, psim] = sim_and_plot(f, w, a, v, p, dt, desc)
    [atest, vtest, ptest] = trapezoidal_int(a(1), v(1,:)', p(1,:)', f, w, dt);
    t = dt * (0:length(p)-1);
    plot_trajectory('traj', ptest)
    plot_error(desc, atest-a, vtest-v, ptest-p, t)
end

function [] = plot_error(desc, ea, ev, ep, t)
    fprintf('%s\n', desc)
    fprintf('velocity error: x1 %.3e, x2 %.3e, abs %.3e\n', max(abs(ev(:,1))), max(abs(ev(:,2))), max(vecnorm(ev,2,2)))
    fprintf('position error: x1 %.3e, x2 %.3e, abs %.3e\n', max(abs(ep(:,1))), max(abs(ep(:,2))), max(vecnorm(ep,2,2)))
    figure
    subplot(3,1,1)
    plot(t, 180/pi*ea')
    ylabel('\alpha [deg]')
    subplot(3,1,2)
    plot(t, ev(:,1), t, ev(:,2))
    ylabel('v [m/s]')
    subplot(3,1,3)
    plot(t, ep(:,1), t, ep(:,2))
    ylabel('x [m]')
    xlabel('t [s]')
    suptitle(desc)
end

function [] = plot_trajectory(name, p)
    figure
    plot(p(:,2), p(:,1))
    ylabel('[m]')
    xlabel('[m]')
    axis('equal')
    title(name)
end

function [a, v, p] = trapezoidal_int(a1, v1, p1, f, w, dt)
    a{1} = a1; v{1} = v1; p{1} = p1;
    for k=2:length(f)
        a{k} = a{k-1} + dt*(w(k-1)+w(k))/2;
        f0 = Rbm(a{k-1})*f(k-1,:)';
        f1 = Rbm(a{k})*f(k,:)';
        v{k} = v{k-1} + (f0 + f1)/2*dt;
        p{k} = p{k-1} + (v{k}+v{k-1})/2*dt;
    end
    a = cell2mat(a)'; v = cell2mat(v)'; p = cell2mat(p)';
end

function R = Rbm(a)
    R = [cos(a) -sin(a); sin(a) cos(a)];
end

function x = GMP(w, dt, beta)
    x = zeros(size(w));
    xk = 0;
    for i = 1:length(w)
        xk = exp(-beta*dt).*xk + w(i,:);
       x(i,:) = xk;
    end
end

function desc = make_description(gyro_rb, gyro_gm, gyro_wn, acc_rb, acc_wn)
    gyro = 'Gyro error:';
    if gyro_rb
        gyro = strcat(gyro, ' Random Bias,');
    end
    if gyro_gm
        gyro = strcat(gyro, ' Gauss-Markov,');
    end
    if gyro_wn
        gyro = strcat(gyro, ' White Noise');
    end
    if ~(gyro_rb || gyro_gm || gyro_wn)
        gyro = strcat(gyro, ' None');
    end
    acc = 'Accel error:';
    if acc_rb
        acc = strcat(acc, ' Random Bias,');
    end
    if acc_wn
        acc = strcat(acc, ' White Noise');
    end
    if ~(acc_rb || acc_wn)
        acc = strcat(acc, ' None');
    end
    desc = sprintf('%s\n%s', gyro, acc);
end
