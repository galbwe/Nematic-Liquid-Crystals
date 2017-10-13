%level set persistence on flat-fielded data
%% import image and apply flat field filter
%image size
image_index = 20009;
[raw_im,ff_im] = FlatFieldFilter(image_index);
fig1 = figure(1)
set(fig1,'position',[0,0,1200,400])
subplot(1,2,1)
imagesc(raw_im)
colormap('gray')
subplot(1,2,2)
imagesc(ff_im)
colormap('gray')
%% compute level set persistence of a window of flat fielded image
[N,M] = size(ff_im);%N is number of pixels in the y direction
yrange = 50:200;%window size
xrange = 50:200;
ff_im2 = ff_im(yrange,xrange);
filtered_persistence = figure(2)
set(filtered_persistence,'position',[0,400,1200,400])
subplot(1,3,1)
imagesc(ff_im2)
title('Image')
colormap('gray')
[H0,H1] = morseFiltration2D(ff_im2);
subplot(1,3,2)
plotpersistencediagram(H0);
title('H_0 persistence')
subplot(1,3,3)
plotpersistencediagram(H1);
title('H_1 persistence')
%% compute persistence of a windowed, unflat-fielded image
[N,M] = size(raw_im);%N is number of pixels in the y direction
yrange = 50:200;%window size
xrange = 50:200;
im2 = raw_im(yrange,xrange);
unfiltered_persistence = figure(3)
set(unfiltered_persistence,'position',[1200,400,1200,400])
subplot(1,3,1)
imagesc(im2)
colormap('gray')
[H0,H1] = morseFiltration2D(im2);
subplot(1,3,2)
plotpersistencediagram(H0);
title('H_0 persistence')
subplot(1,3,3)
plotpersistencediagram(H1);
title('H_1 persistence')
%% generate and compute persistence of a sum of fourier modes
xrange = linspace(0,2*pi,100);
yrange = linspace(0,2*pi,100);
[x,y] = meshgrid(xrange,yrange);
%rows are wavevector direction x, wavevector direction y, amplitude, period
wavedata = [1,1,1,10;-1,1,1,8];
sz = size(wavedata);
number_modes = sz(1,1);
noise = .1*(2*rand(size(x))-1);
z = zeros(size(x));
for i = 1:number_modes
    kx = wavedata(i,1);
    ky = wavedata(i,2);
    ksq = kx^2 + ky^2;
    kx = kx/ksq;
    ky = ky/ksq;
    a = wavedata(i,3);
    p = wavedata(i,4);
    z = z + a*cos(p*(kx*x + ky*y));
end
synthetic_persistence = figure(4)
set(synthetic_persistence,'position',[1200,0,1200,400])
subplot(1,3,1)
imagesc(z)
title('Image')
colormap('gray')
[H0,H1] = morseFiltration2D(z);
subplot(1,3,2)
plotpersistencediagram(H0)
title('H_0 persistence')
subplot(1,3,3)
plotpersistencediagram(H1)
title('H_1 persistence')
%% plot fourier spectrum of the flatfielded data
[N,M] = size(ff_im);
u = ff_im;
uhat = fft2(u);
power = abs(fftshift(uhat));
min_thresh = 2*10^4
max_thresh = 4*10^4
max_peak_height = max(max(power))
[n_0,m_0] = find(power == max_peak_height)
power(power < min_thresh) = 0;
power(power > max_thresh) = max_thresh;
nrange = (n_0-60):(n_0+60);
mrange = (m_0-60):(m_0+60);
nlow = min(nrange);
nhigh = max(nrange);
mlow = min(mrange);
mhigh = max(mrange);
fourier_spectrum = figure(5)
set(fourier_spectrum,'position',[0,800,1200,800])
subplot(2,2,1)
[n,m] = meshgrid(nrange,mrange);
power2 = power(nrange,mrange)
mesh(n,m,power2)
colormap('jet')
xlim([nlow,nhigh])
ylim([mlow,mhigh])
zlim([0,max_thresh])
title('power density spectrum')
subplot(2,2,2)
imagesc(power(nrange,mrange))
%zig mode
subplot(2,2,3)
zig_n_range = 38:44
zig_m_range = 80:93
[zig_n,zig_m] = meshgrid(zig_n_range,zig_m_range);
power3 = power2(zig_n_range,zig_m_range);
mesh(zig_n,zig_m,power3')
colormap('jet')
title('zig modes')
%zag mode
subplot(2,2,4)
zag_n_range = 80:86;
zag_m_range = 78:89;
[zig_n,zig_m] = meshgrid(zag_n_range,zag_m_range);
power4 = power2(zag_n_range,zag_m_range);
%imagesc(power4)
mesh(zig_n,zig_m,power4')
colormap('jet')
title('zag modes')
%% persistence on surface filtered by fourier modes
[N,M] = size(ff_im);
u = ff_im;
uhat = fft2(u);
power = abs(uhat);
min_thresh = 3*10^4;
max_thresh = 4*10^4;
uhat(power < min_thresh) = 0;
uhat(power > max_thresh) = 0;
u2 = real(ifft2(uhat));
filtered_surface = figure(6)
set(filtered_surface,'position',[0,0,1200,1000])
subplot(2,2,1)
imagesc(u2)
colormap('gray')
title('filtered surface')
subplot(2,2,2)
nmin = 300;
nmax = 400;
mmin = 50;
mmax = 130;
nrange = nmin:nmax;
mrange = mmin:mmax;
u3 = u2(nrange,mrange);
imagesc(u3)
title('window for persistence')
[H0,H1] = morseFiltration2D(u3);
subplot(2,2,3)
plotpersistencediagram(H0)
title('H_0 persistence')
subplot(2,2,4)
plotpersistencediagram(H1)
title('H_1 persistence')
%% get the locations of the largest peaks in the fourier spectrum
[N,M] = size(ff_im);
u = ff_im;
uhat = fftshift(fft2(u));
peak_figure = figure(7)
set(peak_figure,'position',[0,0,1600,800])
n = 1;
for num_peaks = 50:50:400
subplot(2,4,n)
n = n+1;
[peaks,i,j] = findpeaks2d(uhat,num_peaks,0,100);
z = zeros(size(uhat));
for k = 1:length(i)
   z(i(k),j(k)) = 1; 
end
imagesc(z);
title(strcat('Number of Peaks = ',num2str(num_peaks)));
end
%% cluster the peaks using kmeans
[N,M] = size(ff_im);
u = ff_im;
uhat = fftshift(fft2(u));
num_peaks = 75;
[peaks,i,j] = findpeaks2d(uhat,num_peaks,0,100);
for t = 1:100
z = zeros(size(uhat));
for k = 1:length(i)
   z(i(k),j(k)) = 1; 
end
figure(8)
imagesc(z);
title(strcat('Number of Peaks = ',num2str(num_peaks)));
data = [i',j'];
num_clusters = 9;
[idX,C] = kmeans(data,num_clusters);
figure(9)
pause(.5)
scatter(C(:,1),C(:,2))
title(strcat('Iteration ',num2str(t)))
end
%% put a window on the zig and zag modes
height = 34;
width = 17;
[N,M] = size(uhat)
C1 = sortrows(C,1);
C2 = sortrows(C1(2:3,:),2);
Czig = floor(C2(1,:))
Czag = floor(C2(2,:))
xzig = floor(Czig(1)-width/2):floor(Czig(1)+width/2-1);
yzig = floor(Czig(2)-height/2):floor(Czig(2)+height/2-1);
xzag = floor(Czag(1)-width/2):floor(Czag(1)+width/2-1);
yzag = floor(Czag(2)-height/2):floor(Czag(2)+height/2-1);
Fzig = zeros(size(uhat));
Fzag = zeros(size(uhat));
Fzig(xzig,yzig) = uhat(xzig,yzig);
Fzag(xzag,yzag) = uhat(xzag,yzag);
Fzag = Fzag + conj(Fzag(N:-1:1,M:-1:1));
Fzig = Fzig + conj(Fzig(N:-1:1,M:-1:1));
Fzag = fftshift(Fzag);
Fzig = fftshift(Fzig);
%convolve with a gaussian
%[y,x] = meshgrid(1:M,1:N);
%a = 10;
%b = 10;
%g = exp(-a*x.^2-b*y.^2);
%Fzig = conv2(g,Fzig,'same');
%Fzag = conv2(g,Fzag,'same');
zig_zag_plots = figure(10)
set(zig_zag_plots,'position',[0,0,1200,800])
colormap('hot')
subplot(2,2,1)
imagesc(abs(fftshift(Fzig)))
title('F_{zig}')
subplot(2,2,2)
imagesc(abs(fftshift(Fzag)))
title('F_{zag}')
subplot(2,2,3)
Izig = real(ifft2(Fzig));
imagesc(Izig)
title('I_{zig}')
subplot(2,2,4)
Izag = real(ifft2(Fzag));
imagesc(Izag)
title('I_{zag}')
figure(11)
imagesc(real(ifft2(Fzig+Fzag)))
colormap('hot')
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
%% compute envelopes of zig and zag modes and plot holes
% shift critical modes to origin implicitly with property of fourier transform
clear
close all
load('mean_power.mat')
M = 640;
N = 480;
T = 30000;
mean_power = fftshift(mean_power2);
figure(13)
imagesc(mean_power)

%load image
image_index = 20009;
[raw_im,I] = FlatFieldFilter(image_index);
F = fft2(I);

%zig modes
M1 = 607;
M2 = 627;
N1 = 16;
N2 = 27;
Wzig = zeros(N,M);
Wzig(N1:N2,M1:M2) = 1;
Fzig = Wzig.*F;
Izig = ifft2(Fzig + conj(Fzig(N:-1:1,M:-1:1)));

%zag modes
M1 = 12;
M2 = 35;
N1 = 17;
N2 = 29;
Wzag = zeros(N,M);
Wzag(N1:N2,M1:M2) = 1;
Fzag = Wzag.*F;
Izag = ifft2(Fzag + conj(Fzag(N:-1:1,M:-1:1)));

%plot zig and zag modes
fig13 = figure(13)
set(fig13,'position',[0,0,1200,400])
subplot(1,3,1)
imagesc(fftshift(real(Izag)));
subplot(1,3,2)
imagesc(fftshift(real(Izig)))
subplot(1,3,3)
imagesc(fftshift(real(Izag + Izig)))
colormap('hot')

%plot envelopes
P_zag_avg = Wzag.*mean_power;
P_zig_avg = Wzig.*mean_power;

%compute critical wavenumbers mc and nc
P_avg = P_zag_avg + P_zig_avg(:,M:-1:1);
figure(15)
imagesc(abs(P_avg))
NP = sum(sum(P_avg(N1:N2,M1:M2)));
[n,m] = meshgrid(M1:M2,N1:N2);
mc = round(sum(sum(m.*P_avg(N1:N2,M1:M2)))/NP);
nc = round(sum(sum(n.*P_avg(N1:N2,M1:M2)))/NP);

%zag envelope - shift (mc,nc) to origin
Azag = 2*real(ifft2(Fzag([nc:N,1:nc-1],[mc:M,1:mc-1])));
%zig envelope - shift (M-mc,nc) to origin
Azig = 2*real(ifft2(Fzig([nc:N,1:nc-1],[M-mc:M,1:M-mc-1])));

fig16 = figure(16)
set(fig16,'position',[0,0,800,400])
subplot(1,2,1)
imagesc(Azag); title('A_{zag}')
subplot(1,2,2)
imagesc(Azig); title('A_{zig}')
colormap('hot')

figure(17)
hold on
plot(1:M,Azag(100,:))
plot(1:M,Izag(100,:))
legend('A_{zag}','I_{zag}')

%find holes (zeros of envelopes)
thresh = 4*1e-2;
Hzag = (abs(Azag) > thresh);
Hzig = (abs(Azig) > thresh);
holes = figure(18);
set(holes,'position',[0,0,800,400]);
subplot(1,2,1)
imagesc(Hzag);
subplot(1,2,2)
imagesc(Hzig);
colormap('gray')
%% plot envelopes and holes over time
clear
close all
load('mean_power.mat')
M = 640;
N = 480;
T = 30000;
imin = 1;
imax = 30000;
thresh = 2*1e-4;
plot_holes = true;

mean_power = fftshift(mean_power2);

%zig window
M1 = 607;
M2 = 627;
N1 = 16;
N2 = 27;
Wzig = zeros(N,M);
Wzig(N1:N2,M1:M2) = 1;

%zag window
M1 = 12;
M2 = 35;
N1 = 17;
N2 = 29;
Wzag = zeros(N,M);
Wzag(N1:N2,M1:M2) = 1;

%compute critical wavenumbers mc and nc
%these will be in the zag window
M1 = 12;
M2 = 35;
N1 = 17;
N2 = 29;
P_zag_avg = Wzag.*mean_power;
P_zig_avg = Wzig.*mean_power;
P_avg = P_zag_avg + P_zig_avg(:,M:-1:1);
NP = sum(sum(P_avg(N1:N2,M1:M2)));
[n,m] = meshgrid(M1:M2,N1:N2);
mc = round(sum(sum(m.*P_avg(N1:N2,M1:M2)))/NP);
nc = round(sum(sum(n.*P_avg(N1:N2,M1:M2)))/NP);

%create a figure to draw on
holes_movie = figure(19)
set(holes_movie,'position',[0,0,1200,800])

for i = imin:imax
    [raw_im,I] = FlatFieldFilter(i);
    F = fft2(I);
    Fzig = Wzig.*F;
    Izig = ifft2(Fzig + conj(Fzig(N:-1:1,M:-1:1)));
    Fzag = Wzag.*F;
    Izag = ifft2(Fzag + conj(Fzag(N:-1:1,M:-1:1)));
    %zag envelope - shift (mc,nc) to origin
    Azag = 2*real(ifft2(Fzag([nc:N,1:nc-1],[mc:M,1:mc-1])));
    %zig envelope - shift (M-mc,nc) to origin
    Azig = 2*real(ifft2(Fzig([nc:N,1:nc-1],[M-mc:M,1:M-mc-1])));
    %holes - zeros of envelopes
    Hzag = (abs(Azag) > thresh);
    Hzig = (abs(Azig) > thresh);
    
    subplot(2,3,1)
    imagesc(fftshift(real(Izig)))
    title('I_{zig}')
    subplot(2,3,2)
    imagesc(fftshift(real(Izag)))
    title('I_{zag}')
    subplot(2,3,3)
    imagesc(fftshift(real(Izig + Izag)))
    title('I_{zig} + I_{zag}')    
    if plot_holes
        subplot(2,3,4)
        imagesc(Hzag);title('H_{zag}');
        subplot(2,3,5)
        imagesc(Hzig);title('H_{zig}');
    else
        subplot(2,3,4)
        imagesc(Azag);title('A_{zag}');
        subplot(2,3,5)
        imagesc(Azig);title('A_{zig}');
    end
    colormap('hot')
    subplot(2,3,6)
    imagesc(I)
    colormap('gray')
    pause(0.5)
end