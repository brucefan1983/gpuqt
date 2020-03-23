clear;font_size=12;close all;
load len_15
load sigma_15;
load len_5
load sigma_5;

figure;
for n=1:5
plot(len_15(1:end,51+n),sigma_15(1:end,51+n),'rs','markersize',4);
hold on;
plot(len_5(1:end,51+n),sigma_5(1:end,51+n),'bo','markersize',4);
p=fminsearch(@(p) norm( p(1)-2/pi*log(len_15(end-30:end,51+n)/p(2)) ...
    - sigma_15(end-30:end,51+n) ), [1,10]);
x=50:500;
plot(x,p(1)-2/pi*log( x/p(2) ), 'k--','linewidth',1);
p=fminsearch(@(p) norm( p(1)+4/pi*log(len_5(end-35:end-20,51+n)/p(2)) ...
    - sigma_5(end-35:end-20,51+n) ), [1,10]);
x=50:500;
plot(x,p(1)+4/pi*log( x/p(2) ), 'k--','linewidth',1);

end
xlabel('$L$ (nm)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
xlim([0,500]);
ylim([0,15]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

text(50,12,'WAL ($\sigma=\sigma_{\rm sc} + 4/\pi \ln(L/l_{\rm sc})$)',...
    'interpreter','latex','color',[0 0 1],'fontsize',12);

text(50,1.5,'WL ($\sigma=\sigma_{\rm sc} - 2/\pi \ln(L/l_{\rm sc})$)',...
    'interpreter','latex','color',[1 0 0],'fontsize',12);

text(300,1.5,'20 meV','fontsize',12, 'interpreter','latex');
text(300,5.2,'100 meV','fontsize',12, 'interpreter','latex');

text(400,5.8,'20 meV','fontsize',12, 'interpreter','latex');
text(400,13,'100 meV','fontsize',12, 'interpreter','latex');

legend('W = 3.77136 eV, \xi = 1.5 a','W = 0.33938 eV, \xi = 5 a');

