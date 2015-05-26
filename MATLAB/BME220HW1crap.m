t = -1:.01:2;
R = 10;
C = 1e-2;
E = .55;
I = zeros(1,301);
I(101:201) = 3e-3;
Vm = (I*R).*(1-exp(-t/(R*C)))+1;
t2 = 1:.01:2
Vm1 = 1+.1*exp(-(t2-1)/(R*C));
%Vm = (I.*R).*(1-exp(t./(R*C)))+E
plot(t,Vm,t2,Vm1)

%%
clear
lam = 2e-3;
tau = 1e-3;
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
t = .000:.0001:.0050;
t = t';
seta1 = Vm(1e-3,t);
seta2 = Vm(1.5e-3,t);
seta3 = Vm(2e-3,t);
% plot(seta1,t,seta2,t,seta3,t),legend('1mm','1.5mm','2mm')
lam = 1e-3;
x = .1e-3:.1e-3:2e-3;
setb1 = Vm(x,.1e-3);
setb2 = Vm(x,.2e-3);
setb3 = Vm(x,.5e-3);
plot(x,setb1,x,setb2,x,setb3),legend('.1ms','.2ms','.5ms')
%%
lam = 1e-3
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
setc1 = Vm(1e-3,t)

lam = 3e-3
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
setc2 = Vm(1e-3,t)

lam = 5e-3
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
setc3 = Vm(1e-3,t)

plot(t,setc1,t,setc2,t,setc3),legend('1','3','5')
%%
tau = 1e-3
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
setc1 = Vm(1e-3,t)

tau = 3e-3
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
setc2 = Vm(1e-3,t)

tau = 5e-3
Vm =@(x,t) (1/(lam.*tau))*sqrt(tau./(pi.*t)).*exp(-((x./lam).^2)./((4*t)/tau)).*exp(-t./tau);
setc3 = Vm(1e-3,t)

plot(t,setc1,t,setc2,t,setc3),legend('1','3','5')