clear;font_size=12;close all;
% data produced by GPUQT
load len_5;
load len_wal;
load len_2;
load len_3;
load len_sqrt3;
load len_15
load sigma_5;
load sigma_wal;
load sigma_3;
load sigma_2;
load sigma_sqrt3;
load sigma_15;

figure;
plot(len_5(1:end-10,51),sigma_5(1:end-10,51),'d','linewidth',1);
hold on;
plot(len_wal(1:end-20,51),sigma_wal(1:end-20,51),'s','linewidth',1);
hold on;
plot(len_3(1:end,51),sigma_3(1:end,51),'o','linewidth',1);
plot(len_2(1:end-20,51),sigma_2(1:end-20,51),'x','linewidth',1);
plot(len_sqrt3(1:end,51),sigma_sqrt3(1:end,51),'+','linewidth',1);
plot(len_15(1:end,51),sigma_15(1:end,51),'v','linewidth',1);
xlabel('$L$ (nm)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
xlim([0,2100]);
ylim([0,10]);

%WAL
x=400:2000;
plot(x,1.42*4+4/pi*log(x/400),'k-','linewidth',2);
plot(x,1.42*4+3/pi*log(x/400),'k-','linewidth',2);

%WL
x=80:2000;
plot(x,4.44-1/pi*log(x/253),'k-','linewidth',2);
[s1,i1]=max(sigma_sqrt3);
plot(x,4-1.8/pi*log(x/143),'k-','linewidth',2);
plot(x,3.33-2/pi*log(x/68),'k-','linewidth',2);

%Sp
x=1:2000;
y=1.42*4*ones(1,2000);
plot(x,y,'k--','linewidth',2);
text(1700,5.2,'$\sigma_{Sp}^{\ast}$','fontsize',12,'interpreter','latex');
text(1000,5,'diffusive');
text(900,6.2,'WAL ($\sigma_{\rm Sp}^{\ast}+3.0/\pi*\ln(L/400 ~\rm nm)$)',...
    'fontsize',12,'interpreter','latex');
text(900,8,'WAL ($\sigma_{\rm Sp}^{\ast}+4.0/\pi*\ln(L/400 ~\rm nm)$)',...
    'fontsize',12,'interpreter','latex');
text(900,3.4,'WL ($4.44-1.0/\pi*\ln(L/253 ~\rm nm)$)',...
        'fontsize',12,'interpreter','latex');
text(900,2.1,'WL ($4.00-1.8/\pi*\ln(L/143~ \rm nm)$)',...
        'fontsize',12,'interpreter','latex');
text(900,0.8,'WL ($3.33-2.0/\pi*\ln(L/68 ~\rm nm)$)',...
        'fontsize',12,'interpreter','latex');

legend('\xi=5a','\xi=4a','\xi=3a','\xi=2a','\xi=1.732a','\xi=1.5a');
