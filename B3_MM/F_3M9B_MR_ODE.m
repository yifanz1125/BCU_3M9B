function dfdt = F_3M9B_MR_ODE(deltacomega)
%% parameters preprocess  
preset = evalin('base','preset');
prefault = evalin('base','prefault');
fault = evalin('base','fault');
postfault = evalin('base','postfault');
Basevalue = evalin('base','Basevalue');
Pm=preset.Pmpu;
E=preset.Epu;
m=preset.m;
d=preset.d;
mT=sum(m,1);
omegab = Basevalue.omegab;

deltac = [deltacomega(1) deltacomega(2) deltacomega(3)];
omega1 = deltacomega(4);
omega2 = deltacomega(5);
omega3 = deltacomega(6);
omega_coi = [omega1 omega2 omega3]*m/mT;


system = evalin('base','system');
switch system
    case "prefault"
        Yred = prefault.Yred;
    case "fault"
        Yred = fault.Yred;
    case "postfault"
        Yred = postfault.Yred;
end

G=real(Yred);
B=imag(Yred);
ngen=size(Yred,1);
Pe=zeros(1,ngen);

%% power calculation
for i=1:ngen
    for j=1:ngen
        delta=deltac(i)-deltac(j);
        Pe(i)=Pe(i)+E(i)*E(j)*(G(i,j)*cos(delta)+B(i,j)*sin(delta));
    end
end


%% differential equation
% deltac
dfdt(1) = omega1 - omega_coi;% deltac1
dfdt(2) = omega2 - omega_coi;% deltac2
dfdt(3) = omega3 - omega_coi;% deltac3

%omega
dfdt(4) = (Pm(1)-Pe(1)-d(1)*(omega1-omegab))/m(1);% omega1
dfdt(5) = (Pm(2)-Pe(2)-d(2)*(omega2-omegab))/m(2);% omega2
dfdt(6) = (Pm(3)-Pe(3)-d(3)*(omega3-omegab))/m(3);% omega3


dfdt = dfdt.';

end