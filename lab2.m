close all;
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

%% Xsens
load('imuData2018/imu5_mtig_c_20181005.mat')
xsens_t = A(:,1);
xsens_Ts = xsens_t(2) - xsens_t(1);
xsens_fIMU = 1/xsens_Ts;
xsens_gyro = A(:,[2, 4]);
xsens_m = mean(xsens_gyro);
xsens_std = std(xsens_gyro);
xsens_N = length(xsens_gyro);

%% Plot data
figure
plot(xsens_t, xsens_gyro(:,1), xsens_t, xsens_gyro(:,2), xsens_t, xsens_gyro(:,3))
legend('gyro x', 'gyro y', 'gyro z')
title('Xsens Gyro')
xlabel('time [s]')
ylabel('[rad/s]')

%% plot PSD
plot_PSD(xsens_gyro, xsens_Ts, 'Xsens')
%% plot autocorrelation
plot_autocorr(xsens_gyro, xsens_Ts, 'Xsens')
%% plot Allan deviation
plot_allandev(xsens_gyro, xsens_Ts, 'Xsens')

% %% plot PSD
% figure
% pwelch(xsens_gyro,[],[],[],1/xsens_Ts,'onesided');
% title('Xsens Power-spectral-density (PSD)')
% legend('x', 'z')
% %% plot autocorrelation
% figure
% g = xsens_gyro - xsens_m;
% [C, LAGS] = xcorr(g, 'unbiased');
% plot(xsens_Ts.*LAGS, C(:,[1,4]))
% xlabel('\tau [s]')
% ylabel('\phi_{xx}(\tau)')
% title('Xsens Autocorrelation (AC)')
% legend('x', 'z')
% %% plot Allan deviation
% figure
% loglog(allandev(xsens_gyro(:,1), '')); hold on; grid on
% %loglog(allandev(xsens_gyro(:,2), ''))
% loglog(allandev(xsens_gyro(:,3), ''))
% title('Allan Deviation')
% xlabel('\tau [s]')
% ylabel('\sigma_y(\tau) [sec]')
% legend('x', 'z')

%% Modeling

% The Xsens IMU does not conatin a Gauss-Markov model.
% autocorrelation is almost perfect Kronecker delta -> beta = inf
% - random bias with mean m
% - white noise with sigma std
% - sinusoid at 41.74 Hz with [-35.49, -54.5] dB/Hz (especially for x axis)
% - perhaps a GM1? No because correlation time would be smaller than Ts

%rng(1);
f_sin = 41.74; % [Hz]
A_sin = [0.0013, 0.00015]; % amplitudes estimated by hand
sinusoid = A_sin.*sin(2*pi*f_sin*xsens_t);
bias = mean(xsens_gyro);
sigma = std(xsens_gyro);
N = length(xsens_gyro);
Ts = xsens_Ts;
model = sigma .* randn(N, 2) + bias + sinusoid;

figure
pwelch(model,[],[],[],1/Ts,'onesided');
title('Xsens Model Power-spectral-density (PSD)')
legend('x', 'z')

figure
g = model - bias;
[C, LAGS] = xcorr(g, 'unbiased');
plot(xsens_Ts.*LAGS, C(:,[1,4]))
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
title('Xsens Model Autocorrelation (AC)')

%% Model 2
beta = [200, 200];
w = sqrt((1-exp(-2*beta*Ts)).*sigma.^2).*randn(N,2);
gm = GMP(w, Ts, beta);
model2 = gm + bias + sinusoid;

% GM models better the slope in the PSD

figure
g = model2 - bias;
[C, LAGS] = xcorr(g, 'unbiased');
plot(xsens_Ts.*LAGS, C(:,[1,4]))
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
title('Xsens Model Autocorrelation (AC)')

figure
pwelch(model2,[],[],[],1/Ts,'onesided');
title('Xsens Model Power-spectral-density (PSD)')
legend('x', 'z')

%% model Allan deviation
figure
loglog(allandev(model2(:,1), '')); hold on; grid on
loglog(allandev(model2(:,2), ''))
title('Allan Deviation')
xlabel('\tau [s]')
ylabel('\sigma_y(\tau) [sec]')
legend('x', 'z')

%% Export data
dlmwrite('02_xsens.txt',xsens_gyro,'precision','%.8f')
dlmwrite('02_xsens_model1.txt',model,'precision','%.8f')
dlmwrite('02_xsens_model2.txt',model2,'precision','%.8f')

%% LN200
[data, fIMU] = readimu('imuData2018/imu2_ln200_b_20181005.imu', 'LN200');
ln200_t = data(:,1) - data(1,1);
ln200_Ts = ln200_t(2) - ln200_t(1);
ln200_fIMU = fIMU;
ln200_gyro = data(:,[2,4]);
ln200_m = mean(ln200_gyro);
ln200_std = std(ln200_gyro);
ln200_N = length(ln200_gyro);
%%
plot_samples(ln200_gyro, ln200_t, 'LN200 Gyro')
%% plot PSD
plot_PSD(ln200_gyro, ln200_Ts, 'LN200')
%% plot autocorrelation
plot_autocorr(ln200_gyro, ln200_Ts, 'LN200')
chi=get(gca, 'Children'); set(gca, 'Children',flipud(chi)) % make x visible
%% plot Allan deviation
plot_allandev(ln200_gyro, ln200_Ts, 'LN200')

%plot_noise_characteristics(ln200_gyro, ln200_t, ln200_Ts, 'LN200 Gyro');
%%
figure
pwelch(ln200_gyro,[],[],[],1/ln200_Ts,'onesided');
title('LN200 Power-spectral-density (PSD)')
%%
figure
plot_autocorr(ln200_gyro, ln200_Ts);
title('LN200 Autocorrelation (AC)')

%% Export data
dlmwrite('02_ln200.txt',ln200_gyro,'precision','%.8f')
dlmwrite('02_ln200_model1.txt',model,'precision','%.8f')
dlmwrite('02_ln200_model2.txt',model2,'precision','%.8f')

%% functions
function x = GMP(w, dt, beta)
    x = zeros(size(w));
    xk = 0;
    for i = 1:length(w)
        xk = exp(-beta*dt).*xk + w(i,:);
       x(i,:) = xk;
    end
end

function [] = plot_samples(samples, t, name)
    figure
    plot(t, samples(:,1), t, samples(:,2))
    legend('x','z')
    title(name)
    xlabel('time [s]')
    ylabel('[rad/s]')
end

function [] = plot_PSD(samples, Ts, name)
    figure
    pwelch(samples,[],[],[],1/Ts,'onesided');
    title(sprintf('%s Power-spectral-density (PSD)',name))
    legend('x', 'z')
end

function [] = plot_autocorr(samples, Ts, name)
    figure
    s = samples - mean(samples);
    [C, LAGS] = xcorr(s, 'unbiased');
    plot(Ts.*LAGS, C(:,[1,4]))
    xlabel('\tau [s]')
    ylabel('\phi_{xx}(\tau)')
    title(sprintf('%s Autocorrelation (AC)',name))
    legend('x', 'z')
end

function [] = plot_allandev(samples, Ts, name)
    figure
    loglog(allandev(samples(:,1), '')); hold on; grid on
    loglog(allandev(samples(:,2), ''))
    title(sprintf('%s Allan Deviation',name))
    xlabel('\tau [s]')
    ylabel('\sigma_y(\tau) [sec]')
    legend('x', 'z')
end
