% Analysis of anisotropic swift hohenberg equation
% u_t = -(1 + d_x^2 + d_y^2)^2u + mu*u + beta*dx^2*u - u^3
% simulate on an LxL grid with periodic boundary conditions
function anisotropicSH()
clc
clear
close all
show_plot = true;
mu = .3;
beta = 2;%Note: there is also a function named Beta, defined at the bottom. 
%Case is used to distinguish between the parameter and the function. Sorry
%if this is excessively confusing.
N = 512;
Nfinal = 100000;
h = 0.1;
ckeep = 100;
L = 30*pi;
filename = 'anisotropicSH';
u = 1e-3*rand(N,N);%random initial condition
ukeep = zeros(int32(Nfinal/ckeep),N,N);
ukeep(1,:,:) = u;
uhat = fft2(u);
tkeep = h*1:ckeep:Nfinal;
ky = pi/L*fftshift(-N/2:N/2-1);
kx = ky;
[kyy,kxx] = meshgrid(kx,ky);
ksq = kxx.^2 + kyy.^2;
Lu = h*(-(1-ksq).^2 + mu - beta*kxx.^2);
Lu2 = .5*Lu;

M = 32; % number of points used in complex contour integral
t = linspace(0,pi,M);
rts = exp(1j*t);
Ikappa = zeros(N,N,M);
Ialpha = zeros(N,N,M);
Ibeta = zeros(N,N,M);
Igamma = zeros(N,N,M);
for m = 1:M
   Ikappa(:,:,m) = Kappa(Lu2 + rts(m));
   Ialpha(:,:,m) = Alpha(Lu + rts(m));
   Ibeta(:,:,m) = Beta(Lu + rts(m));
   Igamma(:,:,m) = Gamma(Lu + rts(m));
end
kappa_u = real(trapz(t,Ikappa,3))/pi;
alpha_u = real(trapz(t,Ialpha,3))/pi;
beta_u = real(trapz(t,Ibeta,3))/pi;
gamma_u = real(trapz(t,Igamma,3))/pi;
eu2 = exp(Lu2);
eu = exp(Lu);

if show_plot
   figs = figure(1)
   set(figs,'position',[0,0,400,700])
end

for n = 1:Nfinal
    
    if show_plot
       subplot(2,1,1)
       imagesc(real(ifft2(uhat)))
       colormap('hot')
       colorbar
       title(strcat('t = ',num2str(n*h)))
       subplot(2,1,2)
       imagesc(fftshift(abs(uhat)))
       pause(.1)
    end
    
    Nu = nonlinear(uhat);
    au = eu2.*uhat + .5*h*kappa_u.*Nu;
    Nau = nonlinear(au);
    bu = eu2.*uhat + .5*h*kappa_u.*Nau;
    Nbu = nonlinear(bu);
    cu = eu2.*au + .5*h*kappa_u.*(2*Nbu - Nu);
    Ncu = nonlinear(cu);
    
    uhat = fft2(real(ifft2(eu.*uhat + h*alpha_u.*Nu + 2*h*beta_u.*(Nau + Nbu) ...
        + h*gamma_u.*Ncu)));
    
    if isnan(norm(uhat))
       'BLOWUP! OH NO!'
       break
    end
    
    if int32(mod(n,ckeep)) == 0
       ukeep(int32(n/ckeep),:,:) = real(ifft2(uhat)); 
       norm(uhat(:),2)
    end
    
end

save(filename,'ukeep','tkeep','beta','mu','N','Nfinal','h','ckeep','L');
end

function k = Kappa(z)
    k = (exp(z) - 1)./z;
end

function a = Alpha(z)
    a = (-4-z + exp(z).*(4-3*z + z.^2))./(z.^3);
end

function b = Beta(z)
    b = (2 + z + exp(z).*(-2+z))./(z.^3);
end

function c = Gamma(z)
    c = (-4 - 3*z - z.^2 + exp(z).*(4-z))./(z.^3);
end

function [uy,ux] = periodic_gradient(u,dx)
    u_ext = zeros(size(u)+4);
    u_ext(3:N+2,3:N+2) = u;
    u_ext(3:N+2,N+3:N+4) = u(:,1:2);
    u_ext(3:N+2,1:2) = u(:,N-1:N);
    u_ext(N+3:N+4,3:N+2) = u(1:2,:);
    u_ext(1:2,3:N+2) = u(N-1:N,:);
    u_ext(N+3:N+4,N+3:N+4) = u(1:2,1:2);
    u_ext(1:2,N+3:N+4) = u(N-1:N,1:2);
    u_ext(N+3:N+4,1:2) = u(1:2,N-1:N);
    u_ext(1:2,1:2) = u(N-1:N,N-1:N);
    [uy_ext,ux_ext] = gradient(u_ext,dx);
    uy = uy_ext(3:N+2,3:N+2);
    ux = ux_ext(3:N+2,3:N+2);
end

function r = rhs(u)
    r = -u.^3;    
end

function nl = nonlinear(uhat)
    u = real(ifft2(uhat));
    u = rhs(u);
    nl = fft2(u);
end

%linear stability analysis
function linearStability(beta,mu)
[ky,kx] = meshgrid(-1.2:0.01:1.2,-1.2:0.01:1.2);
ksq = kx.^2 + ky.^2;
sigma = -(1-ksq).^2 + mu - beta*kx.^2;
figure(1)
imagesc(sigma)
colormap('jet')
colorbar
figure(2)
plot3(kx,ky,sigma)
figure(3)
contour(kx,ky,sigma)
end
