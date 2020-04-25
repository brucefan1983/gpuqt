clear;close all;font_size=13;

% data produced by LSQT
load dos.out;
load msd.out;
dos=(dos+fliplr(dos))/2;
msd=(msd+fliplr(msd))/2;

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
   len(nt,:)=0.142*2*sqrt(msd_ave(nt,:)./dos_ave);
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
for n = 1:9
semilogy(len(:,51+n),sigma_from_msd(:,51+n),'o','color', [1 1 1]*0.05*n,'linewidth',1);
hold on;

p=fminsearch(@(p) norm( p(1)*exp(-len(end-6:end,51+n)/p(2)) - sigma_from_msd(end-6:end,51+n) ),...
    [1,10]);
x=5:30;
semilogy(x,p(1)*exp( -x/p(2) ), 'r--');
end

ylim([0.1,10]);
xlabel('$L$ (nm)', 'fontsize',font_size,'interpreter','latex');
ylabel('$\sigma$ ($e^2/h$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);


text(8,0.12,'0.02 eV','fontsize',font_size,'interpreter','latex');
text(25,0.5,'0.18 eV','fontsize',font_size,'interpreter','latex');

load xi_from_tmm;
load xi_from_sigma;

xi_from_tmm(:,2)=xi_from_tmm(:,2)/2*0.142;
xi_from_sigma(:,2)=xi_from_sigma(:,2)/2*0.142;

axes('Position',[0.45 0.56 0.45 0.33]);
plot(xi_from_tmm(:,1),xi_from_tmm(:,2),'-');
hold on;
plot(xi_from_sigma(28:52,1),xi_from_sigma(28:52,2),'s');
xlim([0,0.19]);
ylim([0,40]);

for n = 1:10
p=fminsearch(@(p) norm( p(1)*exp(-len(end-6:end,51+n)/p(2)) - sigma_from_msd(end-6:end,51+n) ),...
    [1,10]);
plot(n*0.02,p(2),'bv');
end
xlabel('$E$ (eV)','interpreter','latex','fontsize',10);
ylabel('$\xi$ (nm)','interpreter','latex','fontsize',10)
%set(gca,'fontsize',10,'ytick',10.^(0:5));
set(gca,'fontsize',10);
legend('MacKinnon-Kramer','Perturbative estimate of \xi (in 2D)', 'Kubo conductivity scaling');


