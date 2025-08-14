    G=real(Yfull);    
    B=imag(Yfull);
    ngen=size(preset.genno,1);
    nbus=size(Yfull,1);
    E=preset.Epu;
    Id = preset.Id;
    Iq = preset.Iq;    

    Pnet = zeros(nbus-ngen,1);
    Qnet = zeros(nbus-ngen,1);
    delta_net = zeros(nbus-ngen,1);
    V_net = zeros(nbus-ngen,1);
    k=1;
    for i=1:(nbus-ngen)
        delta_net(k)=Results_fsolve(i);
        k=k+1;
    end
    k=1;
    for i=(nbus-ngen)+1:(2*(nbus-ngen))
        V_net(k)=Results_fsolve(i);
        k=k+1;
    end
    clear k


    
    %% P calculation of Bus
    Transform = fault.Transform;
    deltac = delta0;
    pll=pll0;
    for i=1:(nbus-ngen)
        matching_indices = ismember(preset.no_gfl, Transform(i+ngen));
        if any(matching_indices)
            for j=1:ngen
                ddelta=delta_net(i)-deltac(j);
                Pnet(i)=Pnet(i)+E(j)*B(i+ngen,j)*sin(ddelta)+E(j)*G(i+ngen,j)*cos(ddelta);
            end
            for l=1:(nbus-ngen)
                ddelta=delta_net(i)-delta_net(l);
                Pnet(i)=Pnet(i)+V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
            end
            h = find(matching_indices);
            Pnet(i) = Pnet(i)- Id(h)*cos(pll(h)-delta_net(i)) + Iq(h)*sin(pll(h)-delta_net(i));
        else
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
    end
%%
l=1;
i=1;
ddelta=delta_net(i)-delta_net(l);
- Id*cos(pll-delta_net(1)) + Iq*sin(pll-delta_net(1))
-V_net(1)*B(3,3)+ Iq*cos(pll-delta_net(1)) + Id*sin(pll-delta_net(1))

    %% Q calculation of Bus
    for i=1:(nbus-ngen)
        matching_indices = ismember(preset.no_gfl, Transform(i+ngen));
        if any(matching_indices)
            for j=1:ngen
                ddelta=delta_net(i)-deltac(j);
                Qnet(i)=Qnet(i)-E(j)*B(i+ngen,j)*cos(ddelta)+E(j)*G(i+ngen,j)*sin(ddelta);
            end
            for l=1:(nbus-ngen)
                ddelta=delta_net(i)-delta_net(l);
                Qnet(i)=Qnet(i)-V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
            end
            h = find(matching_indices);
            Qnet(i) = Qnet(i) + Iq(h)*cos(pll(h)-delta_net(i)) + Id(h)*sin(pll(h)-delta_net(i));
        else
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
    end
