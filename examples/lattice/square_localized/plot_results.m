clear;font_size=10;

% data produced by GPUQT
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

% now we can define the "true" MSD and the propagating length
for nt=1:Nt
    msd_ave(nt,:)=msd_ave(nt,:)./dos_ave;
end
len=2*sqrt(msd_ave);

% 1D conductance 
g=zeros(Nt,Ne); % g is in unit of 2e^2/h
for ne=1:Ne
   g(:,ne)=sigma_from_msd(:,ne)./len(:,ne)/2;
end
lng=log(g); 

% plot
figure;
subplot(2,2,1);
plot(energy,dos_ave,'k-','linewidth',2);
xlim([-4,4]);
xlabel('Energy ($\gamma$)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/\gamma/a^2$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(a)');

subplot(2,2,2);
surf(energy,t_msd,sigma_from_msd);
xlim([-4,4]);
ylim([0,700]);
zlim([-0.5,5.5]);
xlabel('Energy ($\gamma$)', 'fontsize', font_size,'interpreter','latex');
ylabel('Time ($\hbar/\gamma$)', 'fontsize', font_size,'interpreter','latex');
zlabel('$\sigma$ ($e^2/h$)', 'fontsize', font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
shading interp
title('(b)');

% choose an energy point E and plot VAC(E,t), MSD(E,t), and sigma(E,t)
N=(Ne+1)/2;
p=polyfit(len(5:10,N),lng(5:10,N),1);

subplot(2,2,3);
plot(len(:,N), lng(:,N), 'bo', 'linewidth', 2);
hold on;
x=1:15;
plot(x,polyval(p,x),'r--','linewidth',2);
xlim([1,15]);
ylim([-2,0]);
xlabel('$L$ ($a$)','Fontsize',font_size,'interpreter','latex');
ylabel('$\ln(g)$','Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(c)');

subplot(2,2,4);
plot(t_msd,len(:,N)/pi/2, 'bo', 'linewidth', 2);
hold on;
x=1:800;
plot(x,ones(800,1)*(-0.5/p(1)),'r--','linewidth',2);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('$L/\pi/2$ ($a$)', 'Fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
title('(d)');




