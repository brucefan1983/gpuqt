clear;close all;font_size=12;

%data produced by GPUQT
load 4/dos.out;
load 4/msd.out;
load 4/energy.in;
load 4/time_step.in

% energy points and time steps
Ne=energy(1);
energy=energy(2:Ne+1);
Nt=time_step(1);
time_step=time_step(2:Nt+1);

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
   len(nt,:)=0.142*2*sqrt(msd_ave(nt,:)./dos_ave); % a_CC = 0.142 nm
end

% plot the results
figure;
plot(energy,dos_ave,'linewidth',2);
xlabel('Energy (eV)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/eV/a^2$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

figure;
plot(len(:,51),sigma_from_msd(:,51),'rd','linewidth',1);
xlim([0,2000]);
xlabel('$L$ (nm)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);

len_4=len;
sigma_4=sigma_from_msd;
save('len_4','len_4');
save('sigma_4','sigma_4');




