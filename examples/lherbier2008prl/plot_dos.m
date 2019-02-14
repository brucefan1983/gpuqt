clear; close all; font_size=12;

% load data
data=cell(5,1);
load w=0/dos.out;data{1}=dos;

size(data{1})
load w=1.0/dos.out;data{2}=dos;
load w=1.5/dos.out;data{3}=dos;
load w=2.0/dos.out;data{4}=dos;
load w=2.5/dos.out;data{5}=dos;

% energy points
load w=0/energy.in;
Ne=energy(1);
energy=energy(2:end);

% average over random vectors
for n=1:5
    data{n}=mean(data{n},1);
end

% plot the results for the whole energy range
figure;
plot(energy,data{1},'-','linewidth',2);
hold on;
plot(energy,data{2},'--','linewidth',2);
plot(energy,data{3},':','linewidth',2);
plot(energy,data{4},'-.','linewidth',2);
plot(energy,data{5},'.','linewidth',2);
xlim([-4,4]);
xlabel('Energy ($\gamma$)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/\gamma/A^2$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
legend('ideal','W=1.0\gamma','W=1.5\gamma','W=2.0\gamma','W=2.5\gamma');

% plot the results in the same way as in 
% [Lherbier et al., PRL 100, 036803 (2008)]
energy=energy*2.7; % from gamma to eV
for n=1:5
    data{n}=data{n}/2.7/0.01; % from 1/gamma/A^2 to 1/eV/nm^2
end
figure;
plot(energy,data{1},'-','linewidth',2);
hold on;
plot(energy,data{2},'--','linewidth',2);
plot(energy,data{3},':','linewidth',2);
plot(energy,data{4},'-.','linewidth',2);
plot(energy,data{5},'.','linewidth',2);
xlim([-4,4]);
xlabel('Energy ($eV$)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS ($1/eV/nm^2$)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
legend('ideal','W=1.0\gamma','W=1.5\gamma','W=2.0\gamma','W=2.5\gamma');
