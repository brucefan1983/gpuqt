clear; close all; font_size=10;

% data produced by GPUQT
load dos;
load vac;
load msd;
load sigma_negf;
load sigma_ave;
load sigma_cpgf;


% energy points and time steps
load energy.in;
load time_step.in
Ne=energy(1);
energy=energy(2:Ne+1);
Nt=time_step(1);
time_step=time_step(2:Nt+1);

% average over random vectors
dos_ave=mean(dos,1);
vac_ave=zeros(Nt,Ne);
msd_ave=zeros(Nt,Ne);
Ns=10; % number of independent simulations 
for ns=1:Ns
    index=(ns-1)*Nt+1:ns*Nt;
    vac_ave=vac_ave+vac(index,:);
    msd_ave=msd_ave+msd(index,:);
end
vac_ave=vac_ave/Ns;
msd_ave=msd_ave/Ns;

% conductivity
t_vac=[0;cumsum(time_step(1:end-1))];
t_msd=cumsum(time_step);
sigma_vac=2*pi*cumtrapz(t_vac,vac_ave);
for ne=1:length(energy)
    vac_ave(:,ne)=vac_ave(:,ne)/dos_ave(ne);
end
sigma_msd=zeros(length(time_step),length(energy));
for ne=1:length(energy)
    d_msd=msd_ave(:,ne)-[0;msd_ave(1:end-1,ne)];
    sigma_msd(:,ne)=pi*d_msd./time_step;
    msd_ave(:,ne)=msd_ave(:,ne)/dos_ave(ne);
end

% length and conductance
Ny=50;
L=2*sqrt(msd_ave);

figure;
subplot(2,2,1);
plot(t_vac, vac_ave(:, (Ne+1)/2), 'bo', 'linewidth', 1);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('VAC ($a^2\gamma^2/\hbar^2$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'xtick',0:20:100);
xlim([0,60]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(a)');

subplot(2,2,2);
plot(t_msd, msd_ave(:, (Ne+1)/2), 'rs', 'linewidth', 1);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('MSD ($a^2$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'xtick',0:20:100);
xlim([0,60]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(b)');

subplot(2,2,3);
plot(t_vac(1:30),sigma_vac(1:30,(Ne+1)/2),'bo', 'linewidth', 1);
hold on;
plot(t_msd(1:30), sigma_msd(1:30, (Ne+1)/2), 'rs', 'linewidth', 1);

M=100:100:3000;
eta_scaled=4./M;
eta=eta_scaled*5;
t_cpgf=1./eta;
hold on;
plot(t_cpgf,sigma_ave(:,(Ne+1)/2),'x', 'linewidth', 1);


xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'xtick',0:30:150);
xlim([0,150]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
legend('VAC-KPM','MSD-KPM','KG-CPGF');
% xlabel('hbar/\eta (hbar/\gamma)', 'Fontsize', font_size);
% ylabel('\sigma (e^2/h)', 'Fontsize',font_size);
% set(gca,'xtick',0:30:150);
% xlim([0,150]);
% set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(c)');

subplot(2,2,4);
plot(energy(1:25:end),max(sigma_vac(:,1:25:end)),'bo','linewidth',1);
hold on;
plot(energy(1:25:end),max(sigma_msd(:,1:25:end)),'rs','linewidth',1);
plot(energy(1:25:end),sigma_cpgf(1:25:end),'x','linewidth',1);
plot(sigma_negf(:,1),sigma_negf(:,2),'d','linewidth',1);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('$\sigma_{sc}$ ($e^2/h$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
legend('VAC-KPM','MSD-KPM','KG-CPGF','LB');
ylim([0,70]);
set(gca,'xtick',-5:5);
title('(d)');



