clear;close all;font_size=10;

% data produced by GPUQT
load dos.out;
load vac.out;
load msd.out;

% energy points and time steps
load energy.in;
load time_step.in
Ne=energy(1);
energy=energy(2:end);
Nt=time_step(1);
time_step=time_step(2:end);

% average over random vectors
dos_ave=mean(dos,1);
vac_ave=zeros(Nt,Ne);
msd_ave=zeros(Nt,Ne);
Ns=size(vac,1)/Nt;% number of independent simulations 
for ns=1:Ns
    index=(ns-1)*Nt+1:ns*Nt;
    vac_ave=vac_ave+vac(index,:);
    msd_ave=msd_ave+msd(index,:);
end
vac_ave=vac_ave/Ns;
msd_ave=msd_ave/Ns;

% conductivity from VAC
t_vac = [0; cumsum(time_step(1 : end - 1))];
sigma_from_vac = 2 * pi * cumtrapz(t_vac, vac_ave);

% conductivity from MSD
t_msd=cumsum(time_step)-time_step(1)/2;
sigma_from_msd=zeros(Nt,Ne);
for ne=1:Ne
   sigma_from_msd(:,ne)=pi*(msd_ave(:,ne)-[0;msd_ave(1:end-1,ne)])./time_step;
end

% now we can define the "true" VAC and MSD
for nt=1:Nt
    vac_ave(nt,:)=vac_ave(nt,:)./dos_ave;
    msd_ave(nt,:)=msd_ave(nt,:)./dos_ave;
end

% plot the results
figure;
subplot(3,2,1);
plot(energy,dos_ave,'k-','linewidth',2);
xlim([-6,6]);
xlabel('Energy ($\gamma$)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/\gamma/a^2$)','fontsize',font_size,'interpreter','latex');
title('(a)');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

% choose an energy point E
N=(Ne+1)/2;

subplot(3,2,2);
plot(t_vac,vac_ave(:,N),'rs','linewidth', 2);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('VAC ($a^2\gamma^2/\hbar^2$)', 'Fontsize',font_size,'interpreter','latex');
title('(b)');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

subplot(3,2,3);
plot(t_msd,msd_ave(:,N), 'bo', 'linewidth', 2);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('MSD ($a^2$)', 'Fontsize',font_size,'interpreter','latex');
title('(c)');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

subplot(3,2,4);
plot(t_vac, sigma_from_vac(:,N), 'rs', 'linewidth', 2);
hold on;
plot(t_msd, sigma_from_msd(:,N), 'bo', 'linewidth', 2);
xlabel('Time ($\hbar/\gamma$)','Fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(d)');
legend('from VAC', 'from MSD');

subplot(3,2,5);
surf(energy,t_vac,sigma_from_vac);
xlim([-6,6]);
zlim([0,4]);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('Time ($\hbar/\gamma$)', 'fontsize', font_size,'interpreter','latex');
zlabel('$\sigma$ ($e^2/h$)', 'fontsize', font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(e)');
shading interp

subplot(3,2,6);
plot(energy, max(sigma_from_vac),'r-',energy,max(sigma_from_msd),'b-','linewidth',3);
xlim([-6,6]);
xlabel('Energy ($\gamma$)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma_{sc}$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(f)');
legend('from VAC', 'from MSD')


