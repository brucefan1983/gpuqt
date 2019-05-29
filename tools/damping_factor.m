clear; close all;
M=10000;
m=0:1:M;
a=1/(M+1);

jackson=(1-m*a).*cos(pi*m*a)+a*sin(pi*m*a)/tan(pi*a); % Jackson
lorentz=sinh(4*(1-m/M))/sinh(4);% Lorentz-1
hann=0.5*(1+cos(pi*m*a)); % Hann window

figure;
plot(m,hann,'--','linewidth',2); hold on;
plot(m,jackson,'-.','linewidth',2); hold on;
plot(m,lorentz,'-','linewidth',4);
legend('Hann','Jackson','Lorentz');
ylabel('Damping factor or window function','fontsize',12);
xlabel('m','fontsize',12);
set(gca,'fontsize',12)
set(gca,'ticklength',get(gca,'ticklength')*2)
