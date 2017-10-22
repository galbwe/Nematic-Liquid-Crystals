%% make time indexed arrays with mode values at a single point 
clear
close all
clc
load('mean_power.mat')
mean_power = fftshift(mean_power2);
[N,M] = size(mean_power);
T = 30000;
window1 = [16,30,11,36];
window2 = [454,467,14,36];
cm1 = findCritModes(mean_power,window1);
cm2 = findCritModes(mean_power,window2);
thresh = 0;
l = N/2;
k = M/2;
f1 = zeros(T,1);
f2 = zeros(T,1);
for i = 1:T
    [raw_im,ff_im] = FlatFieldFilter(i);
    I1 = extractSpatialModes(ff_im,window1);
    I2 = extractSpatialModes(ff_im,window2);
    f1(i) = I1(l,k);
    f2(i) = I2(l,k);
end
%% plot time series
load('time_series_10_19_17.mat')
figure(1)
plot(1:T,real(f1))
figure(2)
plot(1:T,real(f2))
%% compute fourier transform of time series data
g1 = fft(f1);
g2 = fft(f2);
figure(3)
plot(-T/2:T/2-1,fftshift(real(g1)))
figure(4)
plot(-T/2:T/2-1,fftshift(real(g2)))
%% compute spatially averaged temporal spectrum
clear
close all
clc
N = 480;
M = 640;
T = 30000;
window1 = [16,30,11,36];
window2 = [454,467,14,36];
%window1 = [1,N,1,M];
%window2 = [1,N,1,M];
Fzag = zeros(window1(2)-window1(1)+1,window1(4)-window1(3)+1,T);
Fzig = zeros(window2(2)-window2(1)+1,window2(4)-window2(3)+1,T);
parfor i = 1:T
    [raw_im,ff_im] = FlatFieldFilter(i);
    ff_im_hat = fft2(ff_im);    
    Fzag(:,:,i) = ff_im_hat(window1(1):window1(2),window1(3):window1(4));
    Fzig(:,:,i) = ff_im_hat(window2(1):window2(2),window2(3):window2(4));
end
save('windowed_time_series','window1','window2','Fzag','Fzig')
%% 
clear
close all
clc
N = 480;
M = 640;
T = 30000;
load('windowed_time_series')
Mw1 = window1(2)-window1(1)+1;
Nw1 = window1(4)-window1(3)+1;
Mw2 = window2(2)-window2(1)+1;
Nw2 = window2(4)-window2(3)+1;
Gzag = fft(Fzag,T,3);
Gzig = fft(Fzig,T,3);
avg_Szag = sum(sum(abs(Gzag).^2))/(Mw1*Nw1);
avg_Szig = sum(sum(abs(Gzig).^2))/(Mw2*Nw2);
temporal_power_spectra = figure(1)
set(temporal_power_spectra,'position',[0,0,1000,400])
subplot(1,2,1)
plot(-T/2:T/2-1,fftshift(avg_Szag(:)),'k')
xlim([-2500,2500])
ylim([0,1e15])
xlabel('w')
ylabel('\langle S_{zag}\rangle')
subplot(1,2,2)
plot(-T/2:T/2-1,fftshift(avg_Szig(:)),'k')
xlim([-2500,2500])
ylim([0,5e14])
xlabel('w')
ylabel('\langle S_{zig}\rangle')
%windows for time spectrum
Wplus = [700,2400];
Wminus = [T-2400,T-700];
SzagPlus = avg_Szag(Wplus(1):Wplus(2));
SzagMinus = avg_Szag(Wminus(2):-1:Wminus(1));
SzigPlus = avg_Szig(Wplus(1):Wplus(2));
SzigMinus = avg_Szig(Wminus(2):-1:Wminus(1));
Sav = SzagPlus + SzagMinus + SzigPlus + SzigMinus;
plot_Sav = figure(2)
plot(Sav(:),'k')
xlabel('w')
ylabel('S_{av}')
%critical Hopf frequency
w = (Wplus(1):Wplus(2))';
wc = round(sum(w.*Sav(:))/sum(Sav(:)));
%filter mask
w = (1:T)';
Mplus = (w > Wplus(1)).*(w < Wplus(2));
Mminus = (w > Wminus(1)).*(w < Wplus(2));
load('critical_modes.mat');
%find temporal envelopes at critical modes
gzagPlus = Gzag(cm1(1)-window1(1)+1,cm1(2)-window1(3)+1,:);
gzagPlus = Mplus.*gzagPlus(:);
gzagMinus = Gzag(cm1(1)-window1(1)+1,cm1(2)-window1(3)+1,:);
gzagMinus = Mminus.*gzagMinus(:);
gzigPlus = Gzig(cm2(1)-window2(1)+1,cm2(2)-window2(3)+1,:);
gzigPlus = Mplus.*gzigPlus(:);
gzigMinus = Gzig(cm2(1)-window2(1)+1,cm2(2)-window2(3)+1,:);
gzigMinus = Mminus.*gzigMinus(:);
fzagPlus = ifft(gzagPlus([wc:T,1:wc-1]));
fzagMinus = ifft(conj(gzagPlus([wc:T,1:wc-1])));
fzigPlus = ifft(gzigPlus([wc:T,1:wc-1]));
fzigMinus = ifft(conj(gzigPlus([wc:T,1:wc-1])));
temporal_envelopes = figure(3)
set(temporal_envelopes,'position',[0,0,800,800])
subplot(4,1,1)
plot(1:T,abs(fzagPlus),'k')
xlabel('t')
ylabel(' |f_{zag}^{+}| ')
subplot(4,1,2)
plot(1:T,abs(fzigPlus),'k')
xlabel('t')
ylabel('| f^{+}_{zig} |')
subplot(4,1,3)
plot(1:T,abs(fzagMinus),'k')
xlabel('t')
ylabel('| f^{-}_{zag} |')
subplot(4,1,4)
plot(1:T,abs(fzigMinus),'k')
xlabel('t')
ylabel('| f^{-}_{zig} |')
%modulated time series
t = (1:T)';
fzagMod = fzagPlus.*exp(2*pi*1j*wc*t/T);
fzagMod = fzagMod + conj(fzagMinus).*exp(-2*pi*1j*wc*t/T);
modulated_time_series = figure(4)
set(modulated_time_series,'position',[0,0,800,800])
t = 20e3:23e3;
subplot(3,1,1)
plot(t,real(fzagPlus(t)),'k')
xlabel('t')
ylabel('Re(f^{+}_{zag})')
subplot(3,1,2)
plot(t,real(fzagMinus(t)),'k')
xlabel('t')
ylabel('Re(f^{-}_{zag})')
subplot(3,1,3)
plot(t,real(fzagMod(t)),'k')
xlabel('t')
ylabel('Re(f^{mod}_{zag})')
%% compute envelopes varying slowly in time and space
Mzag = zeros(size(Gzag));
Mzig = zeros(size(Gzig));
Mzag(:,:,Wplus(1):Wplus(2)) = 1;
Mzig(:,:,Wplus(1):Wplus(2)) = 1;
GzagPlus = Mzag.*Gzag;
GzigPlus = Mzig.*Gzig;
GzagMinus = Mzag.*conj(Gzag);
GzigMinus = Mzig.*conj(Gzig);
FzagPlus = ifft(GzagPlus(:,:,[wc:T,1:wc-1]),T,3);
FzigPlus = ifft(GzigPlus(:,:,[wc:T,1:wc-1]),T,3);
FzagMinus = ifft(GzagMinus(:,:,[wc:T,1:wc-1]),T,3);
FzigMinus = ifft(GzigMinus(:,:,[wc:T,1:wc-1]),T,3);
FzagMinusT = FzagMinus(N:-1:1,M:-1:1,:);
FzigMinusT = FzigMinus(N:-1:1,M:-1:1,:);
nc = cm1(1);
mc = cm1(2);
tmin = 0;
tmax = 30000;
plot_envelopes = figure(5)
set(plot_evelopes,'position',[0,0,800,800])
for i = tmin:tmax
    FzagPlus2 = zeros(N,M);
    FzigPlus2 = zeros(N,M);
    FzagMinus2 = zeros(N,M);
    FzigMinus2 = zeros(N,M);
    FzagPlus2(window1(1):window1(2),window1(3):window1(4)) = FzagPlus(:,:,i);
    FzigPlus2(window2(1):window2(2),window2(3):window2(4)) = FzigPlus(:,:,i);
    FzagMinus2(window1(1):window1(2),window1(3):window1(4)) = FzagMinusT(:,:,i);
    FzigMinus2(window2(1):window2(2),window2(3):window2(4)) = FzigMinusT(:,:,i);
    A1 = ifft2(FzagPlus2([nc:N,1:nc-1],[mc:M,1:mc-1],:));
    A2 = ifft2(FzigPlus2([nc:N,1:nc-1],[mc:M,1:mc-1],:));
    A3 = ifft2(FzagMinus2([nc:N,1:nc-1],[mc:M,1:mc-1],:));
    A4 = ifft2(FzigMinus2([nc:N,1:nc-1],[mc:M,1:mc-1],:));
    subplot(2,2,1)
    imagesc(abs(A1))
    colormap('hot')
    subplot(2,2,2)
    imagesc(abs(A2))
    colormap('hot')
    subplot(2,2,3)
    imagesc(abs(A3))
    colormap('hot')
    subplot(2,2,4)
    imagesc(abs(A4))
    colormap('hot')
    pause(1)
end


