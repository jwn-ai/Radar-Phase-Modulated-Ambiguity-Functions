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

figure,surf(vv,tt,X.^2);
figure,contour(tt, vv, X.^2,100);axis equal;colormap(jet)
X = zeros(2*M-1,Nv);
X0 = zeros(nv,nt,2*M-1);
for p=-(M-1):(M-1)
    i = p+(M-1)+1;
    taup = -p*R-T:dt:T-p*R;
    [ttp, vv] = meshgrid(taup,v);
   
    X0(:,:,i) = (1-abs(ttp)/T).*sinc(T*vv.*(1-abs(ttp)/T));

    figure,surf(vv,ttp,abs(X0(:,:,i)));

    Xv(:,:,i) = (1-abs(p)/M)*abs(diric(2*pi*v*R,M-abs(p)));
    % X(i,:) = Xv(i,:).*X0p;
end

% 
% pi = M-1;
% 
% Xvi = Xv(i,:);
% Xi = X(i,:);
% figure;
% plot(v,Xvi,'--');hold on;
% plot(v,Xi);


% u = single_pulse_ambg(1);
% figure,mesh(u);

function X0 = fX0(tt,vv,T)



end

function [x] = single_pulse_ambg(taup)
eps = 0.000001;
i = 0;
del = 2*taup/150;
for tau = -taup:del:taup
i = i + 1;
j = 0;
fd = linspace(-5/taup,5/taup,151);
val1 = 1. - abs(tau) / taup;
val2 = pi * taup .* (1.0 - abs(tau) / taup) .* fd;
%x(:,i) = abs( val1 .* sin(val2+eps)./(val2+eps));
x(:,i) = abs(val1.*sinc(val2));
T = taup;
[tt,vv] = meshgrid(-taup:del:taup,fd);
tt = abs(tt);
X0 = (1-tt/taup).*sinc(taup*vv.*(1-tt/taup));%

end

end