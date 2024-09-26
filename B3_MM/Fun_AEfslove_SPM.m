% x(2*nbus-ngen) = delta_net(1:nbus-ngen) | V_net(nbus-ngen+1:2*nbus-2*ngen)
function f= Fun_AEfslove_SPM(x,deltac,preset,system)
%system = evalin('base','system');
    prefault = evalin('base','prefault');
    fault = evalin('base','fault');
    postfault = evalin('base','postfault');
    switch system
        case "prefault"
            Yfull = prefault.Yfull_mod;
            Transform = prefault.Transform;
        case "fault1"
            Yfull = fault.Yfull_mod;
            Transform = fault.Transform;
        case "fault2"
            Yfull = fault.Yfull_mod2;
            Transform = fault.Transform2;
        case "postfault"
            Yfull = postfault.Yfull_mod;
            Transform = postfault.Transform;
    end
    G=real(Yfull);    
    B=imag(Yfull);
    ngen=size(preset.genno,1);
    nbus=size(Yfull,1);
    E=preset.Epu;
    

    Pnet = zeros(nbus-ngen,1);
    Qnet = zeros(nbus-ngen,1);
    delta_net = zeros(nbus-ngen,1);
    V_net = zeros(nbus-ngen,1);
    k=1;
    for i=1:(nbus-ngen)
        delta_net(k)=x(i);
        k=k+1;
    end
    k=1;
    for i=(nbus-ngen)+1:(2*(nbus-ngen))
        V_net(k)=x(i);
        k=k+1;
    end
    clear k


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

    for i=1:(nbus-ngen)
        f(i)=Pnet(i);
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

    for i=1:(nbus-ngen)
        f((nbus-ngen)+i)=Qnet(i);
    end

end

