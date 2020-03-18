clear;close all;font_size=12;

% data produced by LSQT
load dos;
load msd;

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
xlim([-5,5]);
xlabel('Energy (eV)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/eV/a^2$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

[sigma_sc,index]=max(sigma_from_msd);

figure;
plot(energy,sigma_sc,'linewidth',2);
xlim([-5,5]);
xlabel('Energy (eV)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma_{max}$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

figure;
semilogx(len(:,1001),sigma_from_msd(:,1001),'rs','linewidth',2);
hold on;
plot(len(:,751),(sigma_from_msd(:,751)+sigma_from_msd(:,1251))/2,'bo','linewidth',2);
x=1:200;
plot(x,sigma_sc(1001)-2/pi*log(x/len(index(1001),1001)),'r-');
plot(x,sigma_sc(751)-2/pi*log(x/len(index(751),751)),'b--');
ylim([0.1,20]);
legend('E = 0 (numerical)','E = \gamma (numerical)','E = 0 (theory)','E = \gamma (theory)');
xlabel('$L$ (a)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

