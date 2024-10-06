% x(2*nbus-ngen-1) = deltac(1: 2) | deltac_net(4:nbus) | V_net(10:9+nbus-ngen)
function dfdt = f_reducedstate2_SPM(x)

postfault = evalin('base','postfault');
preset = evalin('base','preset');

delta2c=x(1);
delta3c=x(2);

Yfull = postfault.Yfull_mod;
ngen=size(preset.genno,1);
nbus=size(Yfull,1);
m=preset.m;
d=preset.d;
Pm=preset.Pmpu;
E=preset.Epu;
mT=sum(m,1);
Pe=zeros(1,ngen);
Pnet = zeros(nbus-ngen,1);
Qnet = zeros(nbus-ngen,1);
Transform = postfault.Transform;
G=real(Yfull);
B=imag(Yfull);

deltac = x(1:2);
delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
deltacc = [delta1c delta2c delta3c];

deltac_net = [x(3) x(4) x(5) x(6) x(7) x(8)];
V_net=[x(9) x(10) x(11) x(12) x(13) x(14)];

%%  power calculation
% P calculation of gen
for i=1:ngen
    for j=1:ngen
        ddelta=deltacc(i)-deltacc(j);
        Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltacc(i)-deltac_net(l);
        Pe(i)=Pe(i)+E(i)*V_net(l)*B(i,l+ngen)*sin(ddelta)+E(i)*V_net(l)*G(i,l+ngen)*cos(ddelta);
    end
end
% P/V calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltacc(j);
        Pnet(i)=Pnet(i)+E(j)*B(i+ngen,j)*sin(ddelta)+E(j)*G(i+ngen,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Pnet(i)=Pnet(i)+V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
               Pnet(i)=Pnet(i)+preset.Sload(h,2)/V_net(i);
        end
    end
end
% Q calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltacc(j);
        Qnet(i)=Qnet(i)-E(j)*B(i+ngen,j)*cos(ddelta)+E(j)*G(i+ngen,j)*sin(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Qnet(i)=Qnet(i)-V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
               Qnet(i)=Qnet(i)+preset.Sload(h,3)/V_net(i);          
        end
    end
end

Pcoi=sum(Pm)-sum(Pe);

dfdt(1) =  (Pm(2)-Pe(2)-Pcoi/sum(m)*m(2))/m(2); 
dfdt(2) =  (Pm(3)-Pe(3)-Pcoi/sum(m)*m(3))/m(3);

% power of load bus
dfdt(3) = -Pnet(1);
dfdt(4) = -Pnet(2);
dfdt(5) = -Pnet(3);
dfdt(6) = -Pnet(4);
dfdt(7) = -Pnet(5);
dfdt(8) = -Pnet(6);
% reactive power of load bus
dfdt(9) = -Qnet(1);
dfdt(10) = -Qnet(2);
dfdt(11) = -Qnet(3);
dfdt(12) = -Qnet(4);
dfdt(13) = -Qnet(5);
dfdt(14) = -Qnet(6);


dfdt = dfdt.';

end