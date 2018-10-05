close all;

rng(1);
N = 200000;
wn = randn(N,3);
rw = cumsum(wn);

%% Gauss-Markov process (GM)

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

for i=1:3
    figure
    hold on;
    plot(wn(:,i))
    plot(rw(:,i))
    plot(gm500(:,i))
    plot(gm2000(:,i))
    legend('WN','RW','GM1, T=500','GM1, T=2000')
    xlabel('sample number []')
    ylabel('signal value []')
end

%% Noise characteristics
% figure
% hold on;
% for i=1:3
%     plot(xcorr(wn(:,i),'unbiased'))
% end
% 
% figure
% hold on;
% for i=1:3
%     plot(xcorr(rw(:,i),'unbiased'))
% end
% 
% figure
% hold on;
% for i=1:3
%     plot(xcorr(gm500(:,i),'unbiased'))
% end

x = gm2000;
figure
title()
subplot(2,1,3)
hold on;
for i=1:3
    plot(xcorr(x(:,i),'unbiased'))
end
subplot(2,1,3)
hold on;
for i=1:3
    plot(pwelch(x(:,i)))
end

%% Autocorrelation
% for i=1:1
%     figure
%     hold on;
%     plot(xcorr(wn(:,i),'unbiased'))
%     plot(xcorr(rw(:,i),'unbiased'))
%     plot(xcorr(gm500(:,i),'unbiased'))
%     plot(xcorr(gm2000(:,i),'unbiased'))
%     legend('WN','RW','GM1, T=500','GM1, T=2000')
%     title('Autocorrelation (AC)')
%     xlabel('sample number []')
%     ylabel('signal value []')
% end

%% Power Spectral Density (PSD)
% for i=1:1
%     figure
%     hold on;
%     plot(pwelch(wn(:,i)))
%     plot(pwelch(rw(:,i)))
%     plot(pwelch(gm500(:,i)))
%     plot(pwelch(gm2000(:,i)))
%     legend('WN','RW','GM1, T=500','GM1, T=2000')
%     title('Power-spectral-density (PSD)')
%     xlabel('sample number []')
%     ylabel('signal value []')
% end

%%
%pwelch(GMP(w, dt, 0.8))

% adev = allandev(gm500(:,1), 'GM T=500');
figure
%plot(adev);
loglog(adev)
grid on
title(sprintf('Allan Deviation: %s', 'GM1, T=500'))
xlabel('\tau [sec]')
ylabel('\sigma_y(\tau) [sec]')
%% functions

function [] = noise_plots(x, desc)
    
end

%     figure
%     
%     plot(AC)
%     title(sprintf('Autocorrelation: %s', name))
%     figure
%     
%     plot(PSD)
%     title(sprintf('Poser Spectral Density: %s', name))

function x = GMP(w, dt, beta)
    x = zeros(size(w));
    xk = 0;
    for i = 1:length(w)
       xk = exp(-beta*dt)*xk + w(i,:);
       x(i,:) = xk;
    end
end
