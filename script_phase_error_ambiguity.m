clear;
close all;

da = .00001;
a = 0:da:1.5;
P = 6;
f = zeros(1,length(a));
A = 1/sqrt(2);
leg_str = {};

%ai = ap(T/2)^p
for p=1:P
    for i=1:length(a)
        ai = a(i);
        fun = @(x) exp(1i*2*pi*ai*x.^p);
        f(i) = .5*abs(integral(fun,-1,1));
    end
    idx = find(abs(f-A) <= da,1);
    a3dB(p) = a(idx);
    figure(1);
    h(p) = plot(a,f,'LineWidth',2);hold on;grid on;
    plot(a3dB(p),A,'g.','MarkerSize',10);
    plot([0,a(end)],[A,A],'--k')
    leg_str = [leg_str, ['p = ' num2str(p)]];
end
legend(h,leg_str);
title('3 dB Phase Error / p^{th} Order Bin Width')
xlabel('\xi_{p}(T/2)^p')
ylabel('Matched Filter Response Magnitude')

a3dB


c1 = 3e8;
c2 = 1500;
v1 = 0:.1:3e6;
v2 = 0:.1:200;
n1 = (c1-v1)./(c1+v1);
n2 = (c2-v2)./(c2+v2);

figure;
semilogx(v1,n1);hold on;
semilogx(v2,n2);

r = 0:.01:.2;
n = (1-r)./(1+r);
figure,plot(r,n);

3e8/4/pi*.1/1000
