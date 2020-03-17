clear;close all;font_size=12;

% data produced by LSQT
load dos.out;
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

% length
len=zeros(Nt,Ne);
for nt=1:Nt
   len(nt,:)=2*sqrt(msd_ave(nt,:)./dos_ave);
end

% plot the results

figure;
surf(energy,t_msd,sigma_from_msd);
xlabel('Energy (eV)', 'fontsize', font_size,'interpreter','latex');
ylabel('Time ($\hbar$/eV)', 'fontsize', font_size,'interpreter','latex');
zlabel('$\sigma$ ($e^2/h$)', 'fontsize', font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
shading interp

figure;
semilogy(energy,dos_ave,'linewidth',2);
xlim([-1,1]);
xlabel('Energy (eV)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/eV/a^2$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

figure;
plot(energy,max(sigma_from_msd),'linewidth',2);
xlim([-1,1]);
xlabel('Energy (eV)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma_{max}$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

figure;
semilogy(len(:,51),sigma_from_msd(:,51),'d','linewidth',2);
hold on;
semilogy(len(:,56),sigma_from_msd(:,56),'o','linewidth',2);
ylim([0.1,10]);
xlabel('L (a)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

