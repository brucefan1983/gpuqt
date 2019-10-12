clear;
close all;

a=1.1;
b=1;
x=0.01:0.01:200;
vac=exp(-x)-0.01*exp(-x/100);
sigma=cumtrapz(x,vac);

figure;
plot(x,sigma,'linewidth',2)
ylim([0,1.1]);
set(gca,'xtick',[],'ytick',[],'fontsize',15);
xlabel('t','fontname','Times New Roman');
ylabel('\sigma(E, t)','fontname','Times New Roman');
text(3,0.3,'I (ballistic)','fontname','Times New Roman','fontsize',15);
text(3,1,'II (diffusive)','fontname','Times New Roman','fontsize',15);
text(90,0.2,'III (localized)','fontname','Times New Roman','fontsize',15);

axes('Position',[.5 .5 .4 .4])
box on

plot(x(1:600),sigma(1:600),'linewidth',2);
hold on;
plot(x(1:95),x(1:95),'--','linewidth',2);
ylim([0,1.2]);
set(gca,'xtick',[],'ytick',[],'fontsize',15);
xlabel('t','fontname','Times New Roman');
ylabel('\sigma(E, t)','fontname','Times New Roman');
text(1,.3,'I (ballistic)','fontname','Times New Roman','fontsize',15);
text(3,1.05,'II (diffusive)','fontname','Times New Roman','fontsize',15);
text(0.1,1.1,'e^2\rho(E)v_F^2(E)t','fontname','Times New Roman','fontsize',12,'color',[1 0 0]);
