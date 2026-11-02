%Phase Modulation Radar Waveforms


clear all;close all;

T = [1e-4 1e-4 1e-4];
B = [100e3 400e3 400e3];
k = -B./T;
dt = 1e-7;
M = 1000;
Ts = M*dt;

for i=2
    t = 0:dt:T(i);
    %t= t-T(i);
s = (exp(j*k(i)*(t).^2));

figure(2)
plot(t,real(s),'LineWidth',2);hold on;grid on;
plot(t,imag(s),'LineWidth',2);
axis([0 T(i) -1 1])
legend('real part','imaginary part')
title(['Up-Chirp Rectangular LFM Pulse T = ' ,num2str(T(i)*1000), ' \mus, B = ', num2str(B(i)/1000) ' kHz'])
xlabel('Time [s]'); ylabel('Amplitude')
z = xcorr(real(s));
figure(1);
plot(-T(i):dt:T(i),z/max(z),'LineWidth',2);hold on;grid on;
title('Compressed Matched Filter Output Pulse')
xlabel('Delay [s]'); ylabel('Amplitude')

end


T = 1e-4;
B = 400e3;
k = B./T;
dt = 1e-6;
Fs = 1/dt;
t = 0:dt:T-dt;
s = real(exp(j*pi*k*t.^2));
s = s/norm(s);
s1 = ones(1,length(t));
s1 = s1/norm(s1);
x = [zeros(1,Fs*4*T) s zeros(1,Fs*2*T)];
x1 = [zeros(1,Fs*4*T) s1 zeros(1,Fs*2*T)];

n = randn(1,length(x));
Ps = .01;
y = x + sqrt(Ps)*n;
y1 = x1 + sqrt(Ps)*n;

z = xcorr(s,y,length(y));
z1 = xcorr(s1,y1,length(y));
figure,plot(y,'LineWidth',2);title('Noisy Signal')
figure;
plot(0:dt:7*T-dt,z(1:length(y)),'LineWidth',2);title('Matched Filter Output')

figure,plot(y1,'LineWidth',2);title('Noisy Signal')
figure;
plot(0:dt:7*T-dt,z1(1:length(y)),'LineWidth',2);title('Matched Filter Output')


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
