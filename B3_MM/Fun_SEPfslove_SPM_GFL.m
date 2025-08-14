% x(2*nbus-ngen) = delta(1: ngen-1) | omegacoi(ngen)| delta_net(ngen+1:nbus) |
% V_net(nbus+1:2*nbus-ngen)| pll(2*nbus-ngen+1: 2*nbus-ngen+ nmachine-ngen) |xint(2*nbus-ngen+ nmachine-ngen+1+nmachine-ngen:2*nbus-ngen+ 2*nmachine-2*ngen)
function f= Fun_SEPfslove_SPM_GFL(x,preset,state,basevalue)
    Y_full=state.Yfull_mod;
    G=real(Y_full);    
    B=imag(Y_full);
    nmachine=size(preset.machineno,1);
    ngen=nmachine-preset.n_gfl;
    nbus=size(Y_full,1);
    m=preset.m;
    d=preset.d;
    Pm=preset.Pmpu;
    E=preset.Epu;
    omegab=basevalue.omegab;
    Id = preset.Id;
    Iq = preset.Iq;
    Kp_pll = preset.Kp_pll;
    Ki_pll = preset.Ki_pll;
    

    delta=zeros(ngen,1); % delta(ngen) is set as 0 as reference
    delta_net = zeros(nbus-ngen,1);
    V_net = zeros(nbus-ngen,1);
    Pe=zeros(ngen,1);
    Pnet = zeros(nbus-ngen,1);
    Qnet = zeros(nbus-ngen,1);
    domegacoi=x(ngen);

    pll = zeros(nmachine-ngen,1); 
    xint = zeros(nmachine-ngen,1); 

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
    k=1;
    for i= (2*nbus-ngen+1):(2*nbus-ngen+nmachine-ngen)
        pll(k)=x(i);
        k=k+1;
    end    
    k=1;
    for i= (2*nbus-ngen+nmachine-ngen+1):(2*nbus-ngen+2*nmachine-2*ngen)
        xint(k)=x(i);
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
        for h=1:preset.n_gfl
            if (preset.no_gfl(h)==state.Transform(i+ngen))
                Pnet(i) = Pnet(i)- V_net(i)*Id(h)*cos(pll(h)-delta_net(i)) + V_net(i)*Iq(h)*sin(pll(h)-delta_net(i));
            end
        end
    end

    for i=1:(nbus-ngen)
        f(ngen+i)=Pnet(i);
    end

    % Q/V calculation of Bus
    for i=1:(nbus-ngen)
        for j=1:ngen
            ddelta=delta_net(i)-delta(j);
            Qnet(i)=Qnet(i)-V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
        end
        for l=1:(nbus-ngen)
            ddelta=delta_net(i)-delta_net(l);
            Qnet(i)=Qnet(i)-V_net(i)*V_net(l)*B(i+ngen,l+ngen)*cos(ddelta)+V_net(i)*V_net(l)*G(i+ngen,l+ngen)*sin(ddelta);
        end
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==state.Transform(i+ngen))
                Qnet(i)=Qnet(i)+preset.Sload(h,3);
            end
        end
        for h=1:preset.n_gfl
            if (preset.no_gfl(h)==state.Transform(i+ngen))
                Qnet(i) = Qnet(i) + V_net(i)*Iq(h)*cos(pll(h)-delta_net(i)) + V_net(i)*Id(h)*sin(pll(h)-delta_net(i));
            end
        end
    end

    for i=1:(nbus-ngen)
        f(nbus+i)=Qnet(i);
    end


    %GFL
    for i=1:(nmachine-ngen)
        pll_bus = preset.no_gfl(i);
        Vq = -1*V_net(pll_bus-ngen)*sin(pll(i)-delta_net(pll_bus-ngen));
        if Ki_pll(i)==0
            f(2*nbus-ngen+i)= Kp_pll(i)*Vq  - domegacoi;
        else
            f(2*nbus-ngen+i)= Kp_pll(i)*Vq + xint(i) - domegacoi;
        end
    end
    for i=1:(nmachine-ngen)
        Vq = -1*V_net(pll_bus-ngen)*sin(pll(i)-delta_net(pll_bus-ngen));
        if Ki_pll(i)==0
            f(2*nbus-2*ngen+nmachine+i)= 0;
        else
            f(2*nbus-2*ngen+nmachine+i)= Ki_pll(i)*Vq;
        end
        
    end


end

