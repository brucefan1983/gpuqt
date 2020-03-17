clear; close all; font_size=12;

% data produced by GPUQT
load dos;
load msd;

% energy points and time steps
load energy.in;
load time_step.in
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

% conductivity
t_msd=cumsum(time_step);
sigma_msd=zeros(length(time_step),length(energy));
for ne=1:length(energy)
    d_msd=msd_ave(:,ne)-[0;msd_ave(1:end-1,ne)];
    sigma_msd(:,ne)=pi*d_msd./time_step;
    msd_ave(:,ne)=msd_ave(:,ne)/dos_ave(ne);
end

% length and conductance
L=2*sqrt(msd_ave);
g_msd=50*sigma_msd./L;
load T_localized;
L_negf=(1:150);
g_negf=T_localized*2;

fminsearch(@(p) norm(p(1)*exp(-L_negf(120:150)/p(2))-g_negf(120:150)),[100,20])

figure;


semilogy(L(1:36, (length(energy)+1)/2), g_msd(1:36, (length(energy)+1)/2),'ro', 'linewidth', 1);
hold on;
semilogy(L_negf(5:5:end),g_negf(5:5:end), 'd','linewidth', 2);
xlabel('$L$ ($a$)', 'Fontsize', font_size,'interpreter','latex');
ylabel('$g$ ($e^2/h$)', 'Fontsize',font_size,'interpreter','latex');
xlim([0,110]);
ylim([0.1,10]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
legend('MSD-KPM','LB');


axes('Position',[0.55 0.55 0.32 0.32])
box on


plot(t_msd, L(:, (length(energy)+1)/2),'ro', 'linewidth', 1);
hold on;
plot(t_msd, 16.5*2*pi*ones(1,length(t_msd)),'--', 'linewidth', 1);
xlabel('Time ($\hbar/\gamma$)', 'Fontsize', 10,'interpreter','latex');
ylabel('$L = 2\sqrt{\Delta X^2(E,t)}$ ($a$)', 'Fontsize',10,'interpreter','latex');
xlim([0,22000]);
ylim([0,130]);
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
text(8000,115,'2\pi\xi(E) from LB','fontsize',10);

