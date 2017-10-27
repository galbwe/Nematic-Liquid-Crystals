%% linear stability
mu = 0.3;
beta = 1;
[ky,kx] = meshgrid(-1.2:0.01:1.2,-1.2:0.01:1.2);
ksq = kx.^2 + ky.^2;
sigma = -(1-ksq).^2 + mu - beta*kx.^2;
figure(1)
imagesc(sigma')
colormap('jet')
colorbar
figure(2)
plot3(kx,ky,sigma)
figure(3)
contour(kx,ky,sigma)
%% plot near-onset SH simulation
filename = '/media/wez/Seagate Backup Plus Drive/nematic-liquid-crystals-8-30-17/SH_simulations/nearOnset-10-24-2017.mat';
load(filename,'Nfinal','ckeep','tkeep');
fig1 = figure(1)
set(fig1,'position',[0,0,1000,400])
iRange = 0:int32(Nfinal/ckeep)-1;
%iRange = 100:100;
for i = iRange
    u = readSHFile(filename,i);
    t = tkeep(i+1);
    uhat = fft2(u);
    subplot(1,2,1)
    imagesc(u')
    colormap('hot')
    colorbar
    title('Anisotropic SH Equation')
    subplot(1,2,2)
    imagesc(fftshift(abs(uhat))) 
    title(['t = ',num2str(t)])
    pause(0.1)
end
%% compute spatially averaged power distribution
clear
close all
clc
filename = '/media/wez/Seagate Backup Plus Drive/nematic-liquid-crystals-8-30-17/SH_simulations/nearOnset-10-24-2017.mat';
load(filename,'Nfinal','ckeep','tkeep','N');
iRange = 50:999;
T = length(iRange);
P_avg = zeros(N,N);
for i = iRange
    u = readSHFile(filename,i);
    uhat = fft2(u);
    P_avg = P_avg + abs(uhat)/T;
end
P_avg = P_avg;
save('../data/SHAvgPower','P_avg')
%% plot spatially averaged power distribution
figure(1)
imagesc(fftshift(P_avg))
colormap('hot')
title('$u^{^}(k_x,k_y)>$')
%% compute critical wavenumbers
clear
close all
clc
filename = '/media/wez/Seagate Backup Plus Drive/nematic-liquid-crystals-8-30-17/SH_simulations/nearOnset-10-24-2017.mat';
load(filename,'Nfinal','ckeep','tkeep','N');
load('../data/SHAvgPower.mat');
%all windowing done in shifted fourier domain
window = [240,280,260,310];%[kymin,kymax,kxmin,kxmax]
shifted_Pavg = fftshift(P_avg);
mask = zeros(N,N);
mask(window(1):window(2),window(3):window(4)) = 1;
masked_Pavg = P_avg.*mask;

%% plot far-from onset SH simulation
filename = '/media/wez/Seagate Backup Plus Drive/nematic-liquid-crystals-8-30-17/SH_simulations/farFromOnset-10-24-2017.mat';
load(filename,'Nfinal','ckeep','tkeep');
fig1 = figure(1)
set(fig1,'position',[0,0,1000,400])
%iRange = 0:int32(Nfinal/ckeep)-1;
iRange = 0:100;
for i = iRange
    u = readSHFile(filename,i);
    t = tkeep(i+1);
    uhat = fft2(u);
    subplot(1,2,1)
    imagesc(u')    
    title('Anisotropic SH Equation')
    colormap('hot')
    colorbar
    subplot(1,2,2)
    imagesc(fftshift(abs(uhat)))
    title(['t = ',num2str(t)])
    pause(0.1)
end