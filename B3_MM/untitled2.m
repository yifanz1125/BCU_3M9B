n=3;
deltac = ep_set(n).xep;
m=preset.m;
delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
deltacc = [delta1c deltac(1) deltac(2)]
theta_net = ep_set(n).delta_net_ep
voltage_net = ep_set(n).voltage_net_ep
f_reducedstate_SPM([ep_set(n).xep;theta_net;voltage_net])
f_reducedstate(ep_set(n).xep)
temp_ini = [theta_net;voltage_net];
system = "postfault";
ccc=Fun_AEfslove_SPM(temp_ini,deltacc,preset,system)
[delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(theta_net+0.01,voltage_net+0.01,deltacc,preset,Basevalue,system,1e5,1e-10)
Results_fsolve1=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc,preset,system),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5))
x_init(ngen)=(postfault.SEP_omegapu-1)*Basevalue.omegab;
for i=1:ngen-1
    x_init(i)=deltacc(i)-deltacc(ngen);
end
for i=1:(nbus-ngen)
    x_init(i+ngen)=theta_net(i)-deltacc(ngen);
end
for i=1:(nbus-ngen)
    x_init(i+nbus)=voltage_net(i);
end
[Results_fsolve2,fval,exitflag,output,jacobian]=fsolve(@(x)Fun_SEPfslove_SPM(x,preset,postfault,Basevalue),x_init,optimset('TolFun',1e-50,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-9));
SEP_omegapu=Results_fsolve2(ngen)/Basevalue.omegab+1
delta_tmp=[Results_fsolve2(1:ngen-1);0];
deltacoi=delta_tmp'*preset.m/sum(preset.m);
SEP_delta=(delta_tmp-deltacoi)
net_delta=Results_fsolve2((ngen+1):nbus)-deltacoi+2*pi
net_voltage=Results_fsolve2((nbus+1):(2*nbus-ngen))


Yfull = postfault.Yfull_mod;
G=real(Yfull);    
B=imag(Yfull);
E=preset.Epu;

Pe=zeros(1,ngen);
Pnet = zeros(nbus-ngen,1);
Qnet = zeros(nbus-ngen,1);
Transform = postfault.Transform;
Pm=preset.Pmpu;
% P calculation of gen
for i=1:ngen
    for j=1:ngen
        ddelta=SEP_delta(i)-SEP_delta(j);
        Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=SEP_delta(i)-net_delta(l);
        Pe(i)=Pe(i)+E(i)*net_voltage(l)*B(i,l+ngen)*sin(ddelta)+E(i)*net_voltage(l)*G(i,l+ngen)*cos(ddelta);
    end
end
% P calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=net_delta(i)-SEP_delta(j);
        Pnet(i)=Pnet(i)+net_voltage(i)*E(j)*B(i+ngen,j)*sin(ddelta)+net_voltage(i)*E(j)*G(i+ngen,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=net_delta(i)-net_delta(l);
        Pnet(i)=Pnet(i)+net_voltage(i)*net_voltage(l)*B(i+ngen,l+ngen)*sin(ddelta)+net_voltage(i)*net_voltage(l)*G(i+ngen,l+ngen)*cos(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
              Pnet(i)=Pnet(i)+preset.Sload(h,2);
        end
    end
end


% Q calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=net_delta(i)-SEP_delta(j);
        Qnet(i)=Qnet(i)-net_voltage(i)*E(j)*B(i+ngen,j)*cos(ddelta)+net_voltage(i)*E(j)*G(i+ngen,j)*sin(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=net_delta(i)-net_delta(l);
        Qnet(i)=Qnet(i)-net_voltage(i)*net_voltage(l)*B(i+ngen,l+ngen)*cos(ddelta)+net_voltage(i)*net_voltage(l)*G(i+ngen,l+ngen)*sin(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
               Qnet(i)=Qnet(i)+preset.Sload(h,3);
        end
    end
end
d=preset.d;
omegab=Basevalue.omegab;
Pm(1)-Pe(1)-d(1)*(SEP_omegapu-1)*omegab
Pm(2)-Pe(2)-d(2)*(SEP_omegapu-1)*omegab
Pm(3)-Pe(3)-d(3)*(SEP_omegapu-1)*omegab