clear;
close all;

t=10.^(-2:0.01:5);
tau_phi=t;
% time evolution approach
VAC=exp(-t);
D=cumtrapz(t,VAC);
MSD=2*cumtrapz(t,D);
% tau_phi regularization
D_phi=tau_phi;
for n=1:length(tau_phi)
    D_phi(n)=trapz(t,exp(-t/tau_phi(n)).*MSD)/2/tau_phi(n)/tau_phi(n);
end

figure;


subplot(2,2,1);
plot(t,D,'r-','linewidth',2)
hold on;
plot(tau_phi,D_phi,'b--','linewidth',2);
ylim([0,1.3]);
xlim([0,50]);
%set(gca,'xtick',[],'ytick',[],'fontsize',12);
xlabel('$t$ or $\tau_{\phi}$','fontname','Times New Roman','interpreter','latex');
ylabel('$D$','interpreter','latex');
text(3,0.3,'I (ballistic)','fontname','Times New Roman','fontsize',12);
text(15,1.1,'II (diffusive)','fontname','Times New Roman','fontsize',12);
legend('D(E,t)',  'D(E,\tau_{\phi})');
title('(a)');

subplot(2,2,2);

plot(t,D,'r-','linewidth',2);
hold on;
plot(tau_phi,D_phi,'b--','linewidth',2);
plot(t,t,'k:','linewidth',2);
ylim([0,0.4]);
xlim([0,0.5]);
%set(gca,'xtick',[],'ytick',[],'fontsize',12);
xlabel('$t$ or $\tau_{\phi}$','fontname','Times New Roman','interpreter','latex');
ylabel('$D$','interpreter','latex');
text(0.1,0.3,'$v_F^2(E)t$','fontsize',12,'color',[0 0 0],'interpreter','latex');
legend('D(E,t)',  'D(E,\tau_{\phi})');
title('(b)');


% time evolution approach
VAC=exp(-t)-0.01*exp(-t/99);
D=cumtrapz(t,VAC);
MSD=2*cumtrapz(t,D);
% tau_phi regularization
D_phi=tau_phi;
for n=1:length(tau_phi)
    D_phi(n)=trapz(t,exp(-t/tau_phi(n)).*MSD)/2/tau_phi(n)/tau_phi(n);
end

subplot(2,2,3);
plot(t,D,'r-','linewidth',2)
hold on;
plot(tau_phi,D_phi,'b--','linewidth',2);
ylim([0,1.1]);
xlim([0,500]);
%set(gca,'xtick',[],'ytick',[],'fontsize',12);
xlabel('$t$ or $\tau_{\phi}$','fontname','Times New Roman','interpreter','latex');
ylabel('$D$','interpreter','latex');
text(3,0.3,'I (ballistic)','fontname','Times New Roman','fontsize',12);
text(10,1,'II (diffusive)','fontname','Times New Roman','fontsize',12);
text(250,0.4,'III (localized)','fontname','Times New Roman','fontsize',12);
legend('D(E,t)',  'D(E,\tau_{\phi})');
title('(c)');

subplot(2,2,4)

plot(t,MSD,'r-','linewidth',2);
hold on;
plot(t,2*tau_phi.*D_phi,'b--','linewidth',2);
xlim([0,5000]);
ylim([0,250]);
%set(gca,'xtick',[],'ytick',[],'fontsize',12);
xlabel('$t$ or $\tau_{\phi}$','fontname','Times New Roman','interpreter','latex');
ylabel('MSD','interpreter','latex');
legend('\DeltaX^2(E,t)', '2\tau_{\phi}D(E,\tau_{\phi})');
title('(d)');
