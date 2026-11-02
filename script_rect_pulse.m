clear; 
T = 1;
tau = -T:0.1:T;
v = -3:.01:3;

% f0 = 1000;
% c = 3e8;
% 
% s = c*f0./v;
% 
% [tt,ss] = meshgrid(tau,v);
% 
% X = (1-abs(tt)/T).*sinc(ss*T.*(1-abs(tt)/T));
% 
% figure,surf(tt,vv,abs(X).^2);
k=1;
for i=1:length(v)
    vi = v(i)
%fun = @(x) ((1-abs(x)/T).*sinc(vi*T.*(1-abs(x)/T))).^2;
fun = @(x) (sin(pi*vi*T.*(1-x/T))/(pi*vi*T)).^2;
        q(k) = integral(fun,0,T);
        k = k+1;
end

figure,plot(v,q)

q0 = 1./(2*(pi*v).^3*T^2).*(v*pi*T-.5*sin(2*pi*v*T));
figure,plot(v,q0)