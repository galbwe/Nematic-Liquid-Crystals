%% creates a signal
clear all;
close all;
clc
L = 10; n = 2048;
t2 = linspace(0,L,n+1);
t = t2(1:n);
k = (2*pi/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k);

S = (3*sin(2*t)+0.5*tanh(0.5*(t-3)) + .2*exp(-(t-4).^2)...
        +1.5*sin(5*t) + 4*cos(3*(t-6).^2))/10+(t/20).^3;
St = fft(S);
%% plot the signal
figure(1)
subplot(2,1,1)
plot(t,S,'k')
set(gca,'Fontsize',14)
xlabel('Time (t)'),ylabel('S(t)')

subplot(2,1,2)
plot(ks,abs(fftshift(St))/max(abs(St)),'k');
axis([-50 50 0 1])
set(gca,'Fontsize',14)
xlabel('frequency(\omega)')
ylabel('FFT(S)')
%% plot Gabor filter with signal
figure(2)
width = [10 1 0.2 2];
for j = 1:4
   g = exp(-width(j)*(t-4).^2);
   subplot(4,1,j)
   plot(t,S,'k'),hold on
   plot(t,g,'k','Linewidth',2)
   set(gca,'FontSize',14)
   ylabel('S(t), g(t)')
end
xlabel('time (t)')
%%
figure(3)
g = exp(-2*(t-4).^2);
Sg=g.*S;
Sgt = fft(Sg);

subplot(3,1,1),plot(t,S,'k'), hold on
plot(t,g,'k','Linewidth',2)
set(gca,'Fontsize',14)
ylabel('S(t), g(t)'), xlabel('time (t)')
subplot(3,1,2), plot(t,Sg,'k')
set(gca,'Fontsize',14)
ylabel('S(t)g(t)'),xlabel('time (t)')
subplot(3,1,3), plot(ks,abs(fftshift(Sgt))/max(abs(Sgt)),'k')
axis([-50,50,0,1])
set(gca,'Fontsize',14)
ylabel('FFT(SG)'),xlabel('frequency (\omega)')