clear;close all;font_size=12;

% data produced by LSQT
load dos.out;

% energy points
load energy.in;
Ne=energy(1);
energy=energy(2:end);

% average over random vectors
dos_ave=mean(dos,1);

% change the units from 1/eV/a^2 to 1/eV/atom
dos_ave=dos_ave*3*sqrt(3)/4;

figure;
plot(energy,dos_ave,'linewidth',2);
xlabel('Energy (eV)','fontsize',font_size,'interpreter','latex');    
ylabel('DOS (1/eV/atom)','fontsize',font_size,'interpreter','latex');
set(gca,'fontsize',font_size,'ticklength',get(gca,'ticklength')*2);
