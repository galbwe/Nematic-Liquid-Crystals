%% compute average of all flat fielded images
M = 640;
N = 480;
T = 30000;
mean_power = zeros(N,M);
for i = 1:T
    [raw_image,u] = FlatFieldFilter(i);
    uhat = fft2(u);
    mean_power = mean_power + abs(uhat)/T;
    i
end
Nrange = N/2-50:N/2+50;
Mrange = M/2-50:M/2+50;
figure(12)
threshold = 10^4;
mean_power2 = fftshift(mean_power);
mean_power2(N/2-10:N/2+10,M/2-10:M/2+10) = 0;
mean_power2(mean_power2 < threshold) = 0;
imagesc(mean_power2(Nrange,Mrange))
save('mean_power.mat','mean_power2')
%% plot mean power and to determine windows and critical wavenumbers
clear
close all
load('mean_power.mat')
M = 640;
N = 480;
mean_power = fftshift(mean_power2);
figure(1)
imagesc(mean_power)
%% compute modes, envelopes and holes
clear
close all
clc
load('mean_power.mat')
mean_power = fftshift(mean_power2);
[N,M] = size(mean_power);
window1 = [16,30,11,36];
window2 = [454,467,14,36];
window3 = [452,466,606,631];
window4 = [15,29,606,629];
%critical modes
%cm1 = [23,23];
%cm2 = [461,23];
%cm3 = [459,619];
%cm4 = [21,619];
cm1 = findCritModes(mean_power,window1);
cm2 = findCritModes(mean_power,window2);
cm3 = findCritModes(mean_power,window3);
cm4 = findCritModes(mean_power,window4);
%threshold for computing holes
thresh = 3e-2;
%times to plot between
imin = 2017;
imax = 2117;
%create figure to draw on
holes_movie = figure(2);
set(holes_movie,'position',[0,0,1600,800]);
for i = imin:imax
    [raw_im,ff_im] = FlatFieldFilter(i);
    [I1,A1,H1] = findEnvAndHoles(ff_im,window1,cm1,thresh);
    [I2,A2,H2] = findEnvAndHoles(ff_im,window2,cm2,thresh);
    [I3,A3,H3] = findEnvAndHoles(ff_im,window3,cm3,thresh);
    [I4,A4,H4] = findEnvAndHoles(ff_im,window4,cm4,thresh);
    %locations of holes
    [X1,Y1] = findOnes(H1);
    [X2,Y2] = findOnes(H2);
    [X3,Y3] = findOnes(H3);
    [X4,Y4] = findOnes(H4);
    %plot envelopes and holes over modes
    subplot(2,4,1)
    imagesc(abs(A1))
    colormap('hot')
    hold on
    scatter(X1,Y1,8,'g','filled')
    xlim([1,M])
    ylim([1,N])
    subplot(2,4,2)
    imagesc(abs(A2))
    colormap('hot')
    hold on
    scatter(X2,Y2,8,'g','filled')
    xlim([1,M])
    ylim([1,N])
    subplot(2,4,3)
    imagesc(abs(A3))
    colormap('hot')
    hold on
    scatter(X3,Y3,8,'g','filled')
    xlim([1,M])
    ylim([1,N])
    subplot(2,4,4)
    imagesc(abs(A4))
    colormap('hot')
    hold on
    scatter(X4,Y4,8,'g','filled')
    xlim([1,M])
    ylim([1,N])
    subplot(2,4,5)
    imagesc(real(I1))
    colormap('hot')
    hold on
    scatter(X1,Y1,8,'g','filled')
    xlim([1,M])
    ylim([1,N])
    subplot(2,4,6)    
    imagesc(real(I2))
    colormap('hot')
    hold on
    scatter(X2,Y2,8,'g','filled')
    subplot(2,4,7)
    imagesc(real(I3))
    colormap('hot')
    hold on
    scatter(X3,Y3,8,'g','filled')
    subplot(2,4,8)
    imagesc(real(I4))
    colormap('hot')
    hold on
    scatter(X4,Y4,8,'g','filled')
    pause(.75)
end