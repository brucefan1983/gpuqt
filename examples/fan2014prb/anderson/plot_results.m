clear;close all;font_size=10;

% data produced by LSQT
load dos.out;
load vac0.out;
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
vac0_ave=mean(vac0,1);
msd_ave=zeros(Nt,Ne);
Ns=size(msd,1)/Nt;% number of independent simulations 
for ns=1:Ns
    index=(ns-1)*Nt+1:ns*Nt;
    msd_ave=msd_ave+msd(index,:);
end
msd_ave=msd_ave/Ns;

% conductivity from MSD
t_msd=cumsum(time_step)-time_step(1)/2;
sigma_from_msd=zeros(Nt,Ne);
for ne=1:Ne
   sigma_from_msd(:,ne)=pi*(msd_ave(:,ne)-[0;msd_ave(1:end-1,ne)])./time_step;
end

% Fermi velocity
velocity=sqrt(vac0_ave./dos_ave);

% mean free path
le=2*max(sigma_from_msd)/2/pi./dos_ave./velocity;

% plot the results

figure;
surf(energy,t_msd,sigma_from_msd);
xlim([-3,3]);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('Time ($\hbar/\gamma$)', 'fontsize', font_size,'interpreter','latex');
zlabel('$\sigma$ ($e^2/h$)', 'fontsize', font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
shading interp

figure;

subplot(2,2,1);
semilogy(energy,dos_ave,'linewidth',2);
xlim([-3,3]);
xlabel('Energy ($\gamma$)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/\gamma/a^2$)','fontsize',font_size,'interpreter','latex');
title('(a)');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

subplot(2,2,2);
plot(energy,max(sigma_from_msd),'linewidth',2);
xlim([-3,3]);
xlabel('Energy ($\gamma$)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma_{sc}$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(b)');

subplot(2,2,3);
plot(energy,velocity,'linewidth',2);
xlim([-3,3]);
xlabel('Energy ($\gamma$)','fontsize',font_size,'interpreter','latex');
ylabel('$v$ ($a\gamma/\hbar$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(c)');

subplot(2,2,4);
semilogy(energy,le,'linewidth',2);
xlim([-3,3]);
xlabel('Energy ($\gamma$)','fontsize',font_size,'interpreter','latex');
ylabel('$l_e$ ($a$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(d)');
