
Yfull = postfault.Yfull_mod;
Transform = postfault.Transform;

G=real(Yfull);    
B=imag(Yfull);
ngen=size(preset.genno,1);
nbus=size(Yfull,1);
E=preset.Epu;
    

Pnet = zeros(nbus-ngen,1);
Qnet = zeros(nbus-ngen,1);
delta_net = zeros(nbus-ngen,1);
V_net = zeros(nbus-ngen,1);
delta_net = delta_net_cri;
V_net=voltage_net_cri;
deltac = deltacc;



    % P calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=delta_net(i)-deltac(j);
        Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=delta_net(i)-delta_net(l);
        Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
    end
    for h=1:size(preset.Sload,1)
       if (preset.Sload(h,1)==Transform(i+ngen))
           if (system == "fault1")||(system == "fault2")
               if Transform(i+ngen)~=fault.faultbus
                    %Pnet(i)=Pnet(i)+preset.Sload(h,2); % during fault
                    %pure impedance
               end
           else
               Pnet(i)=Pnet(i)+preset.Sload(h,2);
           end
       end
    end
end


% Q calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=delta_net(i)-deltac(j);
        Qnet(i)=Qnet(i)-V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=delta_net(i)-delta_net(l);
        Qnet(i)=Qnet(i)-V_net(i)*V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
           if (system == "fault1")||(system == "fault2")
               if Transform(i+ngen)~=fault.faultbus
                    %Qnet(i)=Qnet(i)+preset.Sload(h,3);% during fault
                    %pure impedance
               end
           else
               Qnet(i)=Qnet(i)+preset.Sload(h,3);
           end
           
        end
    end
end

Pnet = zeros(nbus-ngen,1);
Qnet = zeros(nbus-ngen,1);


%%
% delta2c=ep_set(1).xep(1);
% delta3c=ep_set(1).xep(2);
% 
% Yfull = postfault.Yfull_mod;
% ngen=size(preset.genno,1);
% nbus=size(Yfull,1);
% m=preset.m;
% d=preset.d;
% Pm=preset.Pmpu;
% E=preset.Epu;
% mT=sum(m,1);
% Pe=zeros(ngen,1);
% Pnet = zeros(nbus-ngen,1);
% Qnet = zeros(nbus-ngen,1);
% Transform = postfault.Transform;
% G=real(Yfull);
% B=imag(Yfull);
% 
% deltac = ep_set(1).xep;
% delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
% deltacc = [delta1c delta2c delta3c];

deltac = ep_set(1).xep;
delta1c= -preset.m(2:ngen)'*deltac(1:2)/preset.m(1);
deltacc = [delta1c delta2c delta3c];

deltac_net = delta_net_cri;
V_net=voltage_net_cri;

% P calculation of Bus
for i=1:(nbus-ngen)
    for j=1:ngen
        ddelta=deltac_net(i)-deltacc(j);
        Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
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
        ddelta=deltac_net(i)-deltacc(j);
        Qnet(i)=Qnet(i)-V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
    end
    for l=1:(nbus-ngen)
        ddelta=deltac_net(i)-deltac_net(l);
        Qnet(i)=Qnet(i)-V_net(i)*V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
    end
    for h=1:size(preset.Sload,1)
        if (preset.Sload(h,1)==Transform(i+ngen))
               Qnet(i)=Qnet(i)+preset.Sload(h,3);          
        end
    end
end