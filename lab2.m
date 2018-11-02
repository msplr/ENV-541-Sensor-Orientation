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
xsens_allan = plot_allandev(xsens_gyro, xsens_Ts, 'Xsens');

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
[C, LAGS] = xcorr(g, 'coeff');
plot(xsens_Ts.*LAGS, C(:,[1,4]))
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
title('Xsens Model Autocorrelation (AC)')

%% Model 2
beta = [200, 200];
w = sqrt((1-exp(-2*beta*Ts)).*sigma.^2).*randn(N,2);
gm = GMP(w, Ts, beta);
rw = cumsum(1e-6 * randn(N,2));
%model2 = gm  + rw + bias + sinusoid;
model2 = gm;

% GM models better the slope in the PSD

% figure
% g = model2 - bias;
% [C, LAGS] = xcorr(g, 'coeff');
% plot(xsens_Ts.*LAGS, C(:,[1,4]))
% xlabel('\tau [s]')
% ylabel('\phi_{xx}(\tau)')
% title('Xsens Model Autocorrelation (AC)')
% 
% figure
% pwelch(model2,[],[],[],1/Ts,'onesided');
% title('Xsens Model Power-spectral-density (PSD)')
% legend('x', 'z')

%% model Allan deviation
figure
loglog(xsens_allan(:,1)); hold on; grid on
m2_allan = allandev(model2(:,1), '');
loglog(m2_allan)
title('Xsens gyro x compare')
%model2_allan = plot_allandev(m2, Ts, 'Model2');

%% Model 3
f_sin = 41.74; % [Hz]
A_sin = [0.0013, 0.00015]; % amplitudes estimated by hand
sinusoid = A_sin.*sin(2*pi*f_sin*xsens_t);
bias = mean(xsens_gyro);
sigma = std(xsens_gyro);
N = length(xsens_gyro);
Ts = xsens_Ts;

beta = 2.499644e+02;
%w = sqrt((1-exp(-2*beta*Ts)).*sigma.^2).*randn(N,2);
w = sqrt(2.634400e-05) .*randn(N,1);
gm = GMP(w, 1, beta);
wn = sqrt(2.530546e-06) * randn(N, 1);
rw = cumsum(sqrt(2.776097e-12) * randn(N,1));
model3 = wn + rw + gm + bias(1) + sinusoid(:,1);

figure
g = model3 - bias(1);
[C, LAGS] = xcorr(g, 'coeff');
plot(xsens_Ts.*LAGS, C)
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
title('Xsens Model3 Autocorrelation (AC)')

figure
pwelch(model3,[],[],[],1/Ts,'onesided');
title('Xsens Model3 Power-spectral-density (PSD)')
legend('x', 'z')

%% model Allan deviation
model3_allan = plot_allandev(model3, Ts, 'Model3');

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

%% LN200 Model
f_sin = 41.74; % [Hz]
A_sin = [0, 0]; % amplitudes estimated by hand
sinusoid = A_sin.*sin(2*pi*f_sin*ln200_t);
bias = mean(ln200_gyro);
sigma = std(ln200_gyro);
N = length(ln200_gyro);
Ts = ln200_Ts;
model = sigma .* randn(N, 2) + bias + sinusoid;


beta = [200, 200];
w = sqrt((1-exp(-2*beta*Ts)).*sigma.^2).*randn(N,2);
gm = GMP(w, Ts, beta);
model2 = gm + bias + sinusoid;

% GM models better the slope in the PSD

g = model2 - bias;
plot_autocorr(g, ln200_Ts, 'LN200 Model')

plot_PSD(g, ln200_Ts, 'LN200 Model')


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
    [C, LAGS] = xcorr(s, 'coeff');
    plot(Ts.*LAGS, C(:,[1,4]))
    xlabel('\tau [s]')
    ylabel('\phi_{xx}(\tau)')
    title(sprintf('%s Autocorrelation (AC)',name))
    legend('x', 'z')
end

function [adev] = plot_allandev(samples, Ts, name)
    figure
    adev(:,1) = allandev(samples(:,1), '');
    loglog(adev(:,1));
    hold on; grid on
    adev(:,2) = allandev(samples(:,2), '');
    loglog(adev(:,2));
    title(sprintf('%s Allan Deviation',name))
    xlabel('\tau [s]')
    ylabel('\sigma_y(\tau) [sec]')
    legend('x', 'z')
end
