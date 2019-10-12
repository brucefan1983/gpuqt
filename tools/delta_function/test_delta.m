clear; close all;

N=50;
dx=1/25000;
M=10000;

m=0:M;
a=1/(M+1);
jackson=(1-m*a).*cos(pi*m*a)+a*sin(pi*m*a)/tan(pi*a); % Jackson
lorentz=sinh(4*(1-m/M))/sinh(4);% Lorentz-1
hann=0.5*(1+cos(pi*m*a)); % Hann window
eta=0.0001*4;
cpgf=zeros(1,M+1);
for n=0:M
    cpgf(n+1)=(-1)^n*(eta-sqrt(1+eta^2))^n/sqrt(1+eta^2); % CPGF 1
end

T0=zeros(1,M+1);
T0(1)=1;
T0(2)=0;
for n=3:M+1
    T0(n)=-T0(n-2);
end
gamma=[1,2*ones(1,M)];

jackson=jackson.*gamma;
lorentz=lorentz.*gamma;
hann=hann.*gamma;
cpgf=cpgf.*gamma;


kpm_jackson=zeros(1,2*N+1);
kpm_lorentz=zeros(1,2*N+1);
kpm_hann=zeros(1,2*N+1);
kpm_cpgf=zeros(1,2*N+1);

kpm_none=zeros(1,2*N+1);
ftm_none=zeros(1,2*N+1);

ftm_jackson=zeros(1,2*N+1);
ftm_lorentz=zeros(1,2*N+1);
ftm_hann=zeros(1,2*N+1);
ftm_cpgf=zeros(1,2*N+1);


for n=0:2*N
    x=(n-N)*dx;
    Tx=zeros(1,M+1);
    Tx(1)=1;
    Tx(2)=x;
    for m=3:M+1
        Tx(m)=2*x*Tx(m-1)-Tx(m-2);
    end
    
    factor=1/sqrt(1-x^2)/pi;
    kpm_jackson(n+1)=sum(jackson.*T0.*Tx)*factor;
    kpm_lorentz(n+1)=sum(lorentz.*T0.*Tx)*factor;
    kpm_hann(n+1)=sum(hann.*T0.*Tx)*factor;
    kpm_cpgf(n+1)=sum(cpgf.*T0.*Tx)*factor;
    
    ftm_hann(n+1)=sum(cos(pi*(0:M/5)*x).*hann(1:M/5+1))/2;
    
end

x=(-N:N)*dx;

figure;
plot(x,ftm_hann,'*','linewidth',1); hold on;
plot(x,kpm_jackson,'s','linewidth',1); hold on;
plot(x,kpm_lorentz,'o','linewidth',1);
plot(x,kpm_cpgf,'x','linewidth',1);

legend('FTM-Hann','KPM-Jackson','KPM-Lorentz','KPM-CPGF');
ylabel('\delta(x)','fontsize',12);
xlabel('x','fontsize',12);
set(gca,'fontsize',12)
set(gca,'ticklength',get(gca,'ticklength')*2)
