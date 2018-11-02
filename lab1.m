close all;

rng(1);
wn = randn(200000, 3);
rw = cumsum(wn);
% Gauss-Markov process (GM)
dt = 1;
beta = 1/500;
gm500 = GMP(wn, dt, beta);
beta = 1/2000;
gm2000 = GMP(wn, dt, beta);

%% Export data
dlmwrite('01_white_noise.txt',wn,'precision','%.8f')
dlmwrite('01_random_walk.txt',rw,'precision','%.8f')
dlmwrite('01_gm500.txt',gm500,'precision','%.8f')
dlmwrite('01_gm2000.txt',gm2000,'precision','%.8f')

%% Plots
set(groot,'DefaultAxesFontSize',17)
set(groot,'DefaultLineLineWidth',2)

plot_noise_characteristics(wn, 'White Noise')
plot_noise_characteristics(rw, 'Random Walk')
plot_noise_characteristics(gm500, 'Gauss-Markov, T=500')
plot_noise_characteristics(gm2000, 'Gauss-Markov, T=2000')

%%
beta = 1/0.1;
sgm = 1;
dt = 1/100;
wn = sgm*randn(200000, 1);
gm = GMP(wn, dt, beta);
figure
loglog(allandev(gm, ''));
title(sprintf('Guss-Markov, Tc=%f, sgm=%f',1/beta,sgm))
grid on;
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

function [] = plot_noise_characteristics(x, name)
    figure
    subplot(2,2,1);
    plot(x)
    xlabel('sample number []')
    ylabel('signal value []')
    title('Time domain signal')

    subplot(2,2,2);
    tau=-length(x)+1:length(x)-1;
    plot(tau,xcorr(x(:,1), 'coeff')); hold on
    plot(tau,xcorr(x(:,2), 'coeff'))
    plot(tau,xcorr(x(:,3), 'coeff'))
    title('Autocorrelation (AC)')
    xlabel('\tau [s]')
    ylabel('\phi_{xx}(\tau)')

    subplot(2,2,3);
    pwelch(x)
    title('Welch Power-spectral-density (PSD)')

    subplot(2,2,4);
    loglog(allandev(x(:,1), '')); hold on; grid on
    loglog(allandev(x(:,2), ''))
    loglog(allandev(x(:,3), ''))
    title('Allan Deviation')
    xlabel('\tau [s]')
    ylabel('\sigma_y(\tau) [sec]')
    %suptitle(name)
end
