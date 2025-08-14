%% function: calculate potential energy (defined in original system) 
%% Input: thetac_end(n*1)--current point thetac_start--used for path estimation
%% Output: Ep1--generation potential EP2--network potential Ep3-load losses
function [Ep1,Ep2,Ep3,Ep4,Ep5,Ep6,Ep7,Ep8]=Fun_Cal_PotentialEnergy_SPM_GFL(preset,postfault,thetac_end,theta_net_end,voltage_net_end,pll_end,xint_end,omega_coi_end,Basevalue)
    thetac_SEP=postfault.SEP_delta;
    theta_net_SEP=postfault.net_delta;
    voltage_net_SEP=postfault.net_voltage;
    pll_SEP = postfault.pll;
    xint_SEP = postfault.int;
    Yfull = postfault.Yfull_mod;
    G_post=real(Yfull);
    B_post=imag(Yfull);
    ngen=size(preset.genno,1);
    nbus=size(Yfull,1);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    Id = preset.Id;
    Iq = preset.Iq;
    omega_coi_SEP = postfault.SEP_omegapu *Basevalue.omegab;  %coi
    Ep1=0;
    Ep2=0;
    Ep3=0;
    Ep4=0;
    Ep5=0;
    Ep6=0;%pll to net
    Ep7=0;%pll
    Ep8=0;%add


    for i=1:ngen
        Ep1=Ep1+(-1)*(Pm(i)-E(i)^2*G_post(i,i))*(thetac_end(i,1)-thetac_SEP(i,1));
    end
    for i=1:(nbus-ngen)
        for j=1:ngen
            Ep2=Ep2+(-1)*(voltage_net_end(i,1)*E(j)*B_post(i+ngen,j)*cos(theta_net_end(i,1)-thetac_end(j,1))-voltage_net_SEP(i,1)*E(j)*B_post(i+ngen,j)*cos(theta_net_SEP(i,1)-thetac_SEP(j,1)) );
        end
        Ep2=Ep2+(-1)*(voltage_net_end(i,1)^2/2*B_post(i+ngen,i+ngen)-voltage_net_SEP(i,1)^2/2*B_post(i+ngen,i+ngen));
        for h=1:size(preset.Sload,1)
            if (preset.Sload(h,1)==postfault.Transform(i+ngen))
                 Ep5=Ep5+preset.Sload(h,2)*(theta_net_end(i,1)-theta_net_SEP(i,1))+preset.Sload(h,3)*(log(voltage_net_end(i,1))- log(voltage_net_SEP(i,1)));
            end
        end
        matching_indices = ismember(preset.no_gfl, postfault.Transform(i+ngen));
        if any(matching_indices)
            h = find(matching_indices);
            Ep6=Ep6+voltage_net_end(i,1)*Iq(h)*cos(pll_end(h)-thetac_end(i)) + voltage_net_end(i,1)*Id(h)*sin(pll_end(h)-thetac_end(i)) - voltage_net_SEP(i,1)*Iq(h)*cos(pll_SEP(h)-thetac_SEP(i)) - voltage_net_SEP(i,1)*Id(h)*sin(pll_SEP(h)-thetac_SEP(i));
        end
    end
    for i=1:(nbus-ngen-1)
        for l=i+1:(nbus-ngen)
            Ep2=Ep2+(-1)*(voltage_net_end(i,1)*voltage_net_end(l,1)*B_post(i+ngen,l+ngen)*cos(theta_net_end(i,1)-theta_net_end(l,1))-voltage_net_SEP(i,1)*voltage_net_SEP(l,1)*B_post(i+ngen,l+ngen)*cos(theta_net_SEP(i,1)-theta_net_SEP(l,1)) );
        end
    end

    for i=1:(nbus-ngen)
         Ep3=Ep3+G_post(i+ngen,i+ngen)/3*(theta_net_end(i)-theta_net_SEP(i))*( voltage_net_end(i)^2 + voltage_net_end(i)*voltage_net_SEP(i) + voltage_net_SEP(i)^2 );
    end


    % network losses
    if(preset.PathEnergyCal==0)
        Ep4=0;
    elseif(preset.PathEnergyCal==-1)
        Ep4=0;
    else
        n_mid=preset.PathEnergyCal-1;   % numbers of inserted mid-point (n_mid=0--trapezoidal from start to end)  
        dtheta=thetac_end-thetac_SEP;  % dtheta(i)=theta_end(i)-theta_start(i)
        dtheta_net=theta_net_end-theta_net_SEP;
        dvoltage=voltage_net_end-voltage_net_SEP;
        unit_dtheta=dtheta/(n_mid+1);
        unit_dtheta_net=dtheta_net/(n_mid+1);
        unit_dvoltage=dvoltage/(n_mid+1);
        % P lossy of gen
        for i=1:ngen
            for j=1:ngen
                if(i~=j)
                    for m=1:n_mid+1

                        Ep4=Ep4+E(i)*E(j)*G_post(i,j)*0.5*unit_dtheta(i)...,
                            *(cos(thetac_SEP(i)+(m-1)*unit_dtheta(i)-thetac_SEP(j)-(m-1)*unit_dtheta(j))...,
                            +cos(thetac_SEP(i)+m*unit_dtheta(i)-thetac_SEP(j)-m*unit_dtheta(j)));
                    end
                end
            end
            for l=1:(nbus-ngen)
                for m=1:n_mid+1
                    Ep4=Ep4+E(i)*G_post(i,l+ngen)*0.5*unit_dtheta(i)...,
                        *((voltage_net_SEP(l)+(m-1)*unit_dvoltage(l))*cos(thetac_SEP(i)+(m-1)*unit_dtheta(i)-theta_net_SEP(l)-(m-1)*unit_dtheta_net(l))...,
                            +(voltage_net_SEP(l)+m*unit_dvoltage(l))*cos(thetac_SEP(i)+m*unit_dtheta(i)-theta_net_SEP(l)-m*unit_dtheta_net(l)));
                end
            end 
        end
        % P lossy of Bus
        for i=1:(nbus-ngen)
            for j=1:ngen
                for m=1:n_mid+1
                    Ep4=Ep4+E(j)*G_post(i+ngen,j)*0.5*unit_dtheta_net(i)...,
                        *((voltage_net_SEP(i)+(m-1)*unit_dvoltage(i))*cos(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-thetac_SEP(j)-(m-1)*unit_dtheta(j))...,
                        +(voltage_net_SEP(i)+m*unit_dvoltage(i))*cos(theta_net_SEP(i)+m*unit_dtheta_net(i)-thetac_SEP(j)-m*unit_dtheta(j)));
                end
            end
            for l=1:(nbus-ngen)
                if(i~=l)
                    for m=1:n_mid+1
                        Ep4=Ep4+G_post(i+ngen,l+ngen)*0.5*unit_dtheta_net(i)...,
                        *((voltage_net_SEP(i)+(m-1)*unit_dvoltage(i))*(voltage_net_SEP(l)+(m-1)*unit_dvoltage(l))*cos(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-theta_net_SEP(l)-(m-1)*unit_dtheta_net(l))...,
                        +(voltage_net_SEP(i)+(m)*unit_dvoltage(i))*(voltage_net_SEP(l)+(m)*unit_dvoltage(l))*cos(theta_net_SEP(i)+(m)*unit_dtheta_net(i)-theta_net_SEP(l)-(m)*unit_dtheta_net(l)));
                    end
                end
            end
        end
        % Q/V lossy of Bus
        for i=1:(nbus-ngen)
            for j=1:ngen
                for m=1:n_mid+1
                    Ep4 = Ep4+E(j)*G_post(i+ngen,j)*0.5*unit_dvoltage(i)...,
                        *( sin(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-thetac_SEP(j)-(m-1)*unit_dtheta(j))...,
                        + sin(theta_net_SEP(i)+(m)*unit_dtheta_net(i)-thetac_SEP(j)-(m)*unit_dtheta(j)) );
                end
            end
            for l=1:(nbus-ngen)
                if(i~=l)
                    for m=1:n_mid+1
                        Ep4 = Ep4+G_post(i+ngen,l+ngen)*0.5*unit_dvoltage(i)...,
                            *((voltage_net_SEP(l)+(m-1)*unit_dvoltage(l))*sin(theta_net_SEP(i)+(m-1)*unit_dtheta_net(i)-theta_net_SEP(l)-(m-1)*unit_dtheta_net(l)) ...,
                            +(voltage_net_SEP(l)+m*unit_dvoltage(l))*sin(theta_net_SEP(i)+m*unit_dtheta_net(i)-theta_net_SEP(l)-m*unit_dtheta_net(l))  );
                    end
                end
            end
        end
    end
    %pll
    for i=1:preset.n_gfl
        Vq_end = -1*voltage_net_end(preset.no_gfl-ngen)*sin(pll_end(i)-theta_net_end(preset.no_gfl-ngen));
        Vq_SEP = -1*voltage_net_SEP(preset.no_gfl-ngen)*sin(pll_SEP(i)-theta_net_SEP(preset.no_gfl-ngen));
        Ep7 = Ep7 - Vq_end-Vq_SEP; %to be completed


    end


end


