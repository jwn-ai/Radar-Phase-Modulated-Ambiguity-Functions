clear;close all;

% x = -1.5:.01:1.5;
% n = 11;
% p = 2;
% n1 = 5;
% P = 2;
% x2 = x + P;
% s =diric(2*pi*x,n);
% s1 =diric(2*pi*x,n1);
% 
% figure,plot(x,s,'k','LineWidth',2);hold on;grid on;
% xlabel('x=\nuR');
% ylabel(['D_{M}(x) = sin(M\pix)/\pix'])
% title('Dirichlet Kernel for M = 11')
% plot(x,s1,'r','LineWidth',2);
% legend('M=11','M=5')
% 
% 
% pi/(n) %zeros
% 2*pi/(pi*n) % period
% 
% 
% syms v R N;
% 
% f = sin(N*pi*v*R)/sin(pi*v*R)
% 
% simplify(diff(f,1))


f0 = 1e9;
v = 1;
c = 3e8;
v0 = 2.25;
T = .4;

Nv = length(v);
R = 1;
N = 2;
M = 2*N+1;
Td = 2*(N*R + T/2);
dt = .1/5;
dv = .1/5;
v = -v0:dv:v0;
tau0 = -T:dt:T;
tau = -Td:dt:Td;
[tt0,vv0] = meshgrid(tau0,v);
nt = length(tau0);
nv = length(v);
[tt,vv] = meshgrid(tau,v);
X = zeros(size(tt));
for p= -(M-1):(M-1)
    tp = (p*R-T):dt:(p*R+T);
    itp = round(tp/dt + Td/dt + 1);
    iv = round(v/dv + v0/dv + 1);
    Xv = (1-abs(p)/M)*abs(diric(2*pi*vv0*R,M-abs(p)));
X(iv,itp) = Xv.*(1-abs(tt0)/T).*abs(sinc(T*vv0.*(1-abs(tt0)/T)));%
end

figure,surf(tt,vv,X.^1);colormap(jet(16));
xlabel('Delay [s]'); ylabel('Doppler [Hz]');zlabel('|\chi(\tau,\nu)|')
title('Uniform Rectangular Pulse Train Ambiguity for M = 5, T=0.4s, R = 1')

figure,contour(tt, vv, X.^1,20);colormap(jet);grid on;
xlabel('Delay [s]'); ylabel('Doppler [Hz]');
title('Ambiguity Contours')


