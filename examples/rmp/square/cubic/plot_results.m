% figure 9 in the review paper

clear; close all; font_size=13;

% analyzed data produced by GPUQT
load 4000/sigma_e_4000; %1 x Ne
load 3000/sigma_e_3000; %1 x Ne
load 2000/sigma_e_2000; %1 x Ne
load 1000/sigma_e_1000; %1 x Ne
load 3000/sigma_ns_e_3000; %Ns x Ne
Ns=size(sigma_ns_e_3000,1);

% energy points
load 3000/energy.in;
Ne=energy(1);
energy=energy(2:Ne+1);

% main panel
figure
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

% inset
axes('Position',[.35 .24 .4 .4])
plot(energy(1:20:end), sigma_e_1000(1:20:end)-sigma_e_4000(1:20:end), 'd', 'linewidth', 1);
hold on;
plot(energy(1:20:end), sigma_e_2000(1:20:end)-sigma_e_4000(1:20:end), 's', 'linewidth', 1);
plot(energy(1:20:end), sigma_e_3000(1:20:end)-sigma_e_4000(1:20:end), 'o', 'linewidth', 1);
xlim([-7,7]);
ylim([-0.2,0.4]);
xlabel('Energy ($\gamma$)', 'fontsize', font_size/1.3,'interpreter','latex');
ylabel('$\Delta \sigma_{sc}$ ($e^2/ha$)', 'Fontsize',font_size/1.3,'interpreter','latex');
set(gca,'fontsize',font_size/1.3,'ticklength',get(gca,'ticklength')*2);
legend('M=1000 relative to M=4000','M=2000 relative to M=4000','M=3000 relative to M=4000')
