clear; 
T = 1;

v = -3:.01:3;

N = [10 100 1000];  %TB product

leg_str = {};
for j=1:length(N)
Nj = N(j);
k=1;
    for i=1:length(v)
        vi = v(i);
        %fun = @(x) ((1-abs(x)/T).*sinc(vi*T.*(1-abs(x)/T))).^2;

        X2lfm = @(x) ((1-abs(x)).*sinc((Nj*x+vi).*(1-abs(x)))).^2;
        X2qcw = @(x) ((1-abs(x)).*sinc(pi*vi*T.*(1-abs(x)/T))).^2;
        qcw(k) = integral(X2qcw,-T,T);
        qlfm(k) = integral(X2lfm,-T,T);

        k = k+1;
    end
    figure(1);
    plot(v,10*log10(qlfm),'LineWidth',2);hold on;
    leg_str{j} = ['LFM, TB = ' num2str(Nj)];
end

% figure,plot(v,q)
% 
% q0 = 1./(2*(pi*v).^3*T^2).*(v*pi*T-.5*sin(2*pi*v*T));

leg_str{j+1} = 'CW';
figure(1);
plot(v,10*log10(qcw),'LineWidth',2);hold on; grid on;
xlabel('\nuT');
ylabel('Q(\nuT) [dB]');
title('Richman Q-Functions for CW and LFM')
legend(leg_str);