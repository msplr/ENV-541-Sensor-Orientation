close all;

rng(1);
wn = randn(200000, 3);
rw = cumsum(wn);
% Gauss-Markov process (GM)
dt = 1;
beta = 1/2000;
gm500 = GMP(wn, dt, beta);
beta = 1/500;
gm2000 = GMP(wn, dt, beta);

%% Export data
dlmwrite('01_white_noise.txt',wn,'precision','%.8f')
dlmwrite('01_random_walk.txt',rw,'precision','%.8f')
dlmwrite('01_gm500.txt',[rw, gm500, gm2000],'precision','%.8f')
dlmwrite('01_gm500.txt',[rw, gm500, gm2000],'precision','%.8f')

%% Plots
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

i = 1; % examine realization 1
figure; hold on
plot(wn(:,i))
plot(rw(:,i))
plot(gm500(:,i))    
plot(gm2000(:,i))
legend('WN','RW','GM1, T=500','GM1, T=2000')
xlabel('sample number []')
ylabel('signal value []')

%% Autocorrelation
ac_wn = xcorr(wn(:,i), 'unbiased');
ac_rw = xcorr(rw(:,i), 'unbiased');
ac_gm500 = xcorr(gm500(:,i), 'unbiased');
ac_gm2000 = xcorr(gm2000(:,i), 'unbiased');
figure
subplot(2,2,1);
plot(ac_wn)
title('White Noise')
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
subplot(2,2,2);
plot(ac_rw)
title('Random Walk')
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
subplot(2,2,3);
plot(ac_gm500)
title('GM1, T=500')
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
subplot(2,2,4);
plot(ac_gm2000)
title('GM1, T=2000')
xlabel('\tau [s]')
ylabel('\phi_{xx}(\tau)')
suptitle('Autocorrelation (AC)')

%% Power Spectral Density (PSD)
figure; hold on
pwelch([wn(:,i),rw(:,i),gm500(:,i),gm2000(:,i)])
title('Welch Power-spectral-density (PSD)')
legend('WN','RW','GM1, T=500','GM1, T=2000')

%% Allan Deviation
adev_wn = allandev(wn(:,i), '');
adev_rw = allandev(rw(:,i), '');
adev_gm500 = allandev(gm500(:,i), '');
adev_gm2000 = allandev(gm2000(:,i), '');

figure
loglog(adev_wn); hold on; grid on
loglog(adev_rw)
loglog(adev_gm500)
loglog(adev_gm2000)
title('Allan Deviation')
legend('WN','RW','GM1, T=500','GM1, T=2000')
xlabel('\tau [s]')
ylabel('\sigma_y(\tau) [sec]')

%% functions
function x = GMP(w, dt, beta)
    x = zeros(size(w));
    xk = 0;
    for i = 1:length(w)
       xk = exp(-beta*dt)*xk + w(i,:);
       x(i,:) = xk;
    end
end
