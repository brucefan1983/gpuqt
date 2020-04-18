% figure 9 in the review paper

clear; close all; font_size=11;

% analyzed data produced by GPUQT
load sigma_e;
load 3000/sigma_e_3000; %1 x Ne
load 3000/sigma_ns_e_3000; %Ns x Ne
Ns=size(sigma_ns_e_3000,1);

% energy points
load 3000/energy.in;
Ne=energy(1);
energy=energy(2:Ne+1);

% standard errors in the case of M=3000
err3000=zeros(8,1);
for ns=3:10
    err3000(ns-2)=100*mean(std(sigma_ns_e_3000(1:ns,:)))/sqrt(ns)/mean(sigma_e_3000);
end

% relative deviations from M=10 to 400
err10to400=zeros(9,1);
for n=1:9
    err10to400(n)=100*mean(abs(sigma_e(n,:)-sigma_e(end,:)))/mean(sigma_e(end,:));
end

% main panel
figure
subplot(2,2,1:2)

for ns=1:Ns
    plot(energy, sigma_ns_e_3000(ns,:), 'color', [1 1 1]*0.5, 'linewidth', 1);
    hold on;
end
xlim([-7,7]);
ylim([0,15]);
plot(energy, sigma_e_3000, 'r--', 'linewidth', 1);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('$\sigma_{sc}$ ($e^2/ha$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
text(-3,5, '$N = 10^7;~ N_r = 10;~ M = 3000$','interpreter','latex');
title('(a)');

subplot(2,2,3)
plot(3:10,err3000, 'o', 'linewidth', 1);
hold on;
x=1:0.1:11;
plot(x,err3000(4)*sqrt(6)./sqrt(x),'--')
xlim([1,11]);
ylim([0.,1.3]);
xlabel('$N_r$', 'fontsize', font_size,'interpreter','latex');
ylabel('Relative error (%)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
text(2,0.8, '$N_r^{-1/2}$','interpreter','latex');
text(4,1.1, '$N = 10^7;~ M = 3000$','interpreter','latex');
title('(b)');

subplot(2,2,4)
loglog([20:10:50,100:100:400],err10to400(2:end), 'o', 'linewidth', 1);
xlim([10,1000]);
ylim([0.1,20]);
xlabel('$M$', 'fontsize', font_size,'interpreter','latex');
ylabel('Relative error (%)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
text(50,10, '$N = 10^7;~ N_r = 10$','interpreter','latex');
title('(c)');

