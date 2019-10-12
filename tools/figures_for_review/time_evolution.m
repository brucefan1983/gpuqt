clear; close all; font_size=12;

H=0.5;
dt=100;
M=150;

T=zeros(M+1,1);
j=zeros(M+1,1);
U=zeros(M+1,1);

T(1)=1;
T(2)=H;
j(1)=besselj(0,dt);
j(2)=besselj(1,dt);
U(1)=(-1i)^0*j(1)*T(1);
U(2)=U(1)+2*(-1i)^1*j(2)*T(2);

for m=2:M
    j(m+1)=besselj(m,dt);
    T(m+1)=2*H*T(m)-T(m-1);
    U(m+1)=U(m)+2*(-1i)^m*j(m+1)*T(m+1);
end

T=[10:10:90,100:100:900,1000:1000:10000];
count=zeros(length(T),1);
for nt=1:length(T)
    dt=T(nt);
    count(nt)=0;
    for n=1:20000
        count(nt)=count(nt)+1;
        if abs(besselj(n,dt))<1.0e-15
            break;
        end
    end
end
    
figure;

subplot(2,1,1);
plot(1:M+1,real(U),'bs');hold on;
plot(1:M+1,imag(U),'ro');
plot(1:M+1,cos(50)*ones(1,M+1),'b--');
plot(1:M+1,-sin(50)*ones(1,M+1),'r--');
xlabel('N_p','fontsize',font_size);
ylabel('U(\Deltat)','fontsize',font_size);
title('(a)');
text(110,-0.1,'x=0.5; \Deltat=100','fontsize',font_size);
set(gca,'fontsize',font_size);
set(gca,'ticklength',get(gca,'ticklength')*2);
legend('real-numerical','imag-numerical','cos(x\Deltat)','-sin(x\Deltat)');

subplot(2,1,2);
loglog(T,count,'x');hold on;
plot(1:10000,(1:10000)*1.022,'--');
xlabel('\Deltat','fontsize',font_size);
ylabel('N_p','fontsize',font_size);
xlim([10,10000]);
ylim([10,11000]);
title('(b)');
set(gca,'fontsize',font_size);
set(gca,'ticklength',get(gca,'ticklength')*2);
legend('numerical','fit');
