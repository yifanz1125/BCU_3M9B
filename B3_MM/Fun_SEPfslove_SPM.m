% x(2*nbus-ngen) = delta(1: ngen-1) | omegacoi(ngen)| delta_net(ngen+1:nbus) |
% V_net(nbus+1:2*nbus-ngen)
function f= Fun_SEPfslove_SPM(x,preset,state,basevalue)
    Y_full=state.Yfull_mod;
    G=real(Y_full);    
    B=imag(Y_full);
    ngen=size(preset.genno,1);
    nbus=size(Y_full,1);
    m=preset.m;
    d=preset.d;
    Pm=preset.Pmpu;
    E=preset.Epu;
    omegab=basevalue.omegab;
    

    delta=zeros(ngen,1); % delta(ngen) is set as 0 as reference
    Pe=zeros(ngen,1);
    Pnet = zeros(nbus-ngen,1);
    Qnet = zeros(nbus-ngen,1);
    domegacoi=x(ngen);

    for i=1:ngen-1
        delta(i)=x(i);
    end
    k=1;
    for i=ngen+1:nbus
        delta_net(k)=x(i);
        k=k+1;
    end
    k=1;
    for i=nbus+1:(2*nbus-ngen)
        V_net(k)=x(i);
        k=k+1;
    end
    clear k

    % Pe calculation
    for i=1:ngen
        for j=1:ngen
            ddelta=delta(i)-delta(j);
            Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta(i)-delta_net(l);
            Pe(i)=Pe(i)+E(i)*V_net(l)*B(i,l+ngen)*sin(ddelta)+E(i)*V_net(l)*G(i,l+ngen)*cos(ddelta);
        end
    end

    PCOI=sum(Pm-Pe);   
    
    for i=1:ngen-1
        f(i)=Pm(i)-Pe(i)-m(i)/sum(m)*PCOI+m(i)/sum(m)*sum(d)*domegacoi-d(i)*domegacoi;
    end
    f(ngen)= sum(Pm-Pe)-sum(d)*domegacoi;    % the nth equ

    % P calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-delta(j);
            Pnet(i)=Pnet(i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Pnet(i)=Pnet(i)+V_net(i)*V_net(l)*B(i+ngen,l+ngen)*sin(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*cos(ddelta);
        end
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==state.Transform(i+ngen))
                Pnet(i)=Pnet(i)+preset.Sload(h,2);
            end
        end
    end

    for i=1:(nbus-ngen)
        f(ngen+i)=Pnet(i);
    end

    % Q calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-delta(j);
            Qnet(i)=Qnet(i)-E(j)*B(i+ngen,j)*cos(ddelta)+E(j)*G(i+ngen,j)*sin(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Qnet(i)=Qnet(i)-V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
        end
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==state.Transform(i+ngen))
                Qnet(i)=Qnet(i)+preset.Sload(h,3)/V_net(i);
            end
        end
    end

    for i=1:(nbus-ngen)
        f(nbus+i)=Qnet(i);
    end

end

