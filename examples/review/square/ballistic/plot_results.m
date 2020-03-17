clear; close all; font_size=12;

% data produced by GPUQT
load dos.out;
load vac.out;
load msd.out;

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
Ns=size(vac,1)/Nt % number of independent simulations 
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

% conductance
L=2*sqrt(msd_ave);
g_msd=2*sigma_msd./L;
g_vac=zeros(length(time_step),length(energy));
for ne=1:length(energy)
    g_vac(:,ne)=dos_ave(ne)*sqrt(vac_ave(:,ne));
end
g_vac=g_vac*2*pi;
load T;
E_negf=-4.1:0.01:4.1;
g_negf=T*2;

figure;
plot(energy,mean(g_vac),'linewidth',2);
hold on;
plot(E_negf,g_negf,'r--','linewidth',2);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('Conductance $g(E)$ ($e^2/h$)', 'fontsize', font_size,'interpreter','latex');
ylim([0, 10]);
set(gca,'xtick',-5:5);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
legend('VAC-KPM (MSD-KPM)','LB');

font_size=10;
figure;

subplot(2,2,1);
plot(t_vac, vac_ave(:, (length(energy)+1)/2), 'bo', 'linewidth', 1);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('VAC ($a^2\gamma^2/\hbar^2$)', 'Fontsize',font_size,'interpreter','latex');
ylim([0,4]);
set(gca,'xtick',0:2:10);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(a)', 'Fontsize',font_size);

subplot(2,2,2);
plot(cumsum(time_step), msd_ave(:, (length(energy)+1)/2), 'bo', 'linewidth', 1);
hold on;
plot(0:0.1:10,3*(0:0.1:10).^2,'k-','linewidth',2);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('MSD ($a^2$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'xtick',0:2:10);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(b)', 'Fontsize',font_size);

subplot(2,2,3);
plot(energy, dos_ave, 'r-', 'linewidth', 2);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('DOS ($1/\gamma/a^2$)', 'fontsize', font_size,'interpreter','latex');
xlim([-4, 4]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(c)', 'Fontsize',font_size);

vac_ave(:,1:250)=0;
vac_ave(:,end-249:end)=0;
subplot(2,2,4);
plot(energy, sqrt(vac_ave(1,:)), 'r-', 'linewidth', 2);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('$v$ ($a\gamma/\hbar$)', 'fontsize', font_size,'interpreter','latex');
xlim([-4, 4]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(d)', 'Fontsize',font_size);

