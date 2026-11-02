close all; clear all;

%ambiguity function
N = 50;  %TB product
dx = .001;
dy = .005;
x = -1:dx:1;  %x = tau/T
y = -3:dy:3;  %y = v*T;
[xx,yy] = meshgrid(x,y);
axx = abs(xx);
X = @(xx,yy) (1-abs(xx)).*sinc((N*xx+yy).*(1-abs(xx)));

figure;
subplot(2,2,1);
surf(xx,yy,abs(X(xx,yy)),'EdgeColor','none');colormap(jet);
xlabel('\tau/T'); ylabel('\nuT');
zlabel('|\chi(\tau,\nu)|');
title(['LFM Ambiguity for TB = ' num2str(N)])

subplot(2,2,2);
contour(xx, yy, abs(X(xx,yy)));colormap(jet);grid on;
xlabel('\tau/T'); ylabel('\nuT');
title('LFM Ambiguity Contours')

subplot(2,2,3);
plot(y,abs(X(0,y)),'LineWidth',2);grid on;
xlabel('\nuT');ylabel('|\chi(0,\nu)|');
title('LFM Ambiguity Doppler Profile')

subplot(2,2,4);
plot(x,abs(X(x,0)),'LineWidth',2);grid on;
xlabel('\tau/T');ylabel('|\chi(\tau,0)|');
title('LFM Ambiguity Delay Profile')


%%%%Spectrum

%z = f/B;

N = [10 100 1000];
z = 0:.001:1;
legend_str = {};
for i=1:length(N)

t1 = sqrt(2*N(i))*(.5-z);
t2 = sqrt(2*N(i))*(.5+z);



C1 = fresnelc(t1);
C2 = fresnelc(t2);
S1 = fresnels(t1);
S2 = fresnels(t2);

M = (C1+C2).^2 + (S1+S2).^2;
M = M/max(M(:));

figure(3)
plot(z,M,'LineWidth',1);grid on;hold on;
legend_str{i} = ['TB = ' num2str(N(i))];
end
xlabel('f/B');ylabel('|M(f)|');
title('Frequency Spectrum of Rectangular Chirp')
legend(legend_str);
