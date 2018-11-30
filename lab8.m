close all
w = pi/100; % [rad]
r = 25; % [m]
t_end = round(2*pi/w); % [s]
dt = 1; % [s] Kalman filter update time
t_KF = (dt:dt:t_end)'; % [s] timestamp of KF updates

% reference trajectory
a_ref = @(t) t*w + pi/2;
v_ref = @(t) (w*r*[cos(a_ref(t)) sin(a_ref(t))]);
p_ref = @(t) (r*[sin(a_ref(t)) -cos(a_ref(t))]);

%% Measurements
dt_gps = 1; % [s]
t_gps = (dt_gps:dt_gps:t_end)'; % [s] timestamp of GPS measurements
std_gps = [0.5 0.5]; % [m] GPS standard deviation
z_gps = p_ref(t_gps) + std_gps.*randn(length(t_gps),2);

fprintf('Empirical std dev GPS:%.4f\n', sqrt(sum(var(z_gps-p_ref(t_gps)))))
%% Kalman Filter
model = 'a=const';

if model == 'a=const'
    F = diag([1,1,1,1],2);
    G = [zeros(4,2); eye(2)];
    H = [eye(2), zeros(2,4)];
    W = diag([0.01^2 0.01^2]); % [(m/s^3/sqrt(Hz)]^2
    f0 = [-r*w^2 0];
    x0 = [p_ref(0) v_ref(0) f0]';
    P0 = diag([10^2 10^2 0.1^2 0.1^2 0.1^2 0.1^2]);
elseif model == 'v=const'
    F = diag([1,1],2);
    G = [zeros(2); eye(2)];
    H = [eye(2), zeros(2)];
    W = diag([0.05^2, 0.05^2]); % [(m/s^2/sqrt(Hz)]^2
    x0 = [p_ref(0) v_ref(0)]';
    P0 = diag([10^2 10^2 0.1^2 0.1^2]);
else
    fprintf('ERROR: unknown model\n')
    return
end
R = diag(std_gps.^2);
[Phi, Qk] = KF_Phi_and_Qk(F, G, W, dt);

% vectors tracking evolution
nx = length(Phi);
x_KF_pred = zeros(length(t_KF),nx);
x_KF = zeros(length(t_KF),nx);
sigma_KFp = zeros(length(t_KF),1);
innovation = zeros(length(t_gps),2);

x = x0;
P = P0;
i_gps = 1;
for i = 1:length(t_KF)
    [xp, Pp] = KF_predict(x,P,Phi,Qk);
    if t_KF(i) >= t_gps(i_gps)
        z = z_gps(i_gps,:)';
        [x, P] = KF_correct(xp,Pp,z,H,R);
        innovation(i_gps,:) = z - H*xp;
        i_gps = i_gps + 1;
    else
        x = xp;
        P = Pp;
    end
    x_KF_pred(i,:) = xp';
    x_KF(i,:) = x';
    sigma_KFp(i) = sqrt(P(1,1) + P(2,2)); % KF-predicted positioning quality
end

fprintf('Empirical std dev KF filter:%.4f\n', sqrt(sum(var(x_KF(:,1:2)-p_ref(t_KF)))))
fprintf('Final std dev KF predicted:%.4f\n', sigma_KFp(end))

%% Plots
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

plot_traj(p_ref(t_KF), x_KF, x_KF_pred, z_gps)

figure
plot(t_KF, sigma_KFp)
title('KF-predicted positioning quality')
xlabel('time [s]'); ylabel('\sigma^{KFp}_{xy} [m]')

t_stable_P = 25; % [s]
plot_signal_and_hist(innovation, t_gps, 'All innovations', 'innovation y [m]', 'y_{north}', 'y_{east}')
plot_signal_and_hist(innovation(t_gps>t_stable_P,:), t_gps(t_gps>t_stable_P), ...
    'Innovations stable gain', 'innovation y [m]', 'y_{north}', 'y_{east}')
plot_signal_and_hist(x_KF(:,3:4)-v_ref(t_KF), t_KF, 'Velocity error', '[m/s]', 'v_{err} (north)', 'v_{err} (east)')

%% save plots
%saveplot(figure(1), '08/traj.pdf')
saveplot(figure(2), '08/sigma_KFp.pdf')
saveplot(figure(3), '08/innovation_all.pdf', 2, 1)
saveplot(figure(4), '08/innovation_stable.pdf', 2, 1)
saveplot(figure(5), '08/velocity_error.pdf', 2, 1)

%% functions
function [x, P] = KF_predict(x,P,Phi,Q)
    x = Phi*x;
    P = Phi*P*Phi' + Q;
end

function [x, P] = KF_correct(x,P,z,H,R)
    K = P*H'*(H*P*H' + R)^-1;
    x = x + K*(z - H*x);
    P = (eye(size(P)) - K*H)*P;
end

function [Phi, Qk] = KF_Phi_and_Qk(F, G, W, dt)
    n = length(F);
    A = [-F, G*W*G'; ...
              zeros(n,n), F'];
    B = expm(dt*A);
    Phi = B(n+1:end,n+1:end)';
    Qk = Phi*B(1:n,n+1:end);
end

function [] = plot_signal_and_hist(innovation, t, name, yl, l1, l2)
    figure
    ax1 = subplot(1,2,1);
    plot(t, innovation(:,1)); hold on
    plot(t, innovation(:,2))
    legend(l1, l2)
    xlabel('time [s]'); ylabel(yl)
    set(gca,'XLim',[t(1) t(end)])
    title(name)
    ax2 = subplot(1,2,2);
    edges=linspace(min(min(innovation)), max(max(innovation)), 10);
    histogram(innovation(:,1), edges,'Orientation', 'horizontal'); hold on
    histogram(innovation(:,2), edges,'Orientation', 'horizontal')
    xlabel('counts')
    % same y axis
    set(ax2,'ytick',[])
    ylim = get([ax1, ax2], {'YLim'}); ylim = cat(2, ylim{:});
    set([ax1, ax2], 'Ylim', [min(ylim), max(ylim)])
    % adjust subplot ratio
    p1 = get(ax1, 'Position');
    p2 = get(ax2, 'Position');
    gap = p2(1)-p1(1)-p1(3);
    dx = 0.5*p1(3);
    set(ax1, 'Position', [p1(1:2), p1(3)+dx+0.5*gap, p1(4)])
    set(ax2, 'Position', [p2(1)+dx, p2(2), p1(3)-dx, p1(4)])
end

function [] = plot_traj(p, x_KF, x_KF_pred, z_gps)
    figure;
    plot(p(:,2), p(:,1)); hold on;
    plot(x_KF(:,2), x_KF(:,1))
    plot(x_KF_pred(:,2), x_KF_pred(:,1))
    plot(z_gps(:,2), z_gps(:,1), 'x')
    title('Trajectory'); xlabel('x2 [m]'); ylabel('x1 [m]')
    legend('true position', 'corrected position', 'predicted position', 'GPS measurement')
    axis equal
end

function saveplot(h,outfilename,xscale,yscale)
    if nargin < 4
        yscale = 1;
    end
    if nargin < 3
        xscale = 1;
    end
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [xscale*pos(3) yscale*pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 xscale*pos(3) yscale*pos(4)]);
    print(outfilename, '-dpdf');
end