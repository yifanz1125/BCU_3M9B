%% function: check the err of SEP result
%% input: omega_SEP (in SRF frame, rad/s)
function [Perr,flag_err]=Fun_SEPcheck(State,preset,delta_SEP,omega_SEP)
    Yred=State.Yred;
    ngen=size(Yred,1);
    B=imag(Yred);
    G=real(Yred);
    E=preset.Epu;
    d=preset.d;
    m=preset.m;
    Pm=preset.Pmpu;
    flag_err=0;
    Pe_tmp=zeros(ngen,1);
    for i=1:ngen
        for j=1:ngen
            ddelta=delta_SEP(i)-delta_SEP(j);
            Pe_tmp(i)=Pe_tmp(i)+E(i)*E(j)*(B(i,j)*sin(ddelta)+G(i,j)*cos(ddelta));
        end
    end
    Pcoi_tmp=sum(Pm)-sum(Pe_tmp)-sum(d)*omega_SEP;
    Err_PCUEP=Pm-Pe_tmp-Pcoi_tmp/sum(m)*m-d*omega_SEP;
    Perr=Err_PCUEP;
    if(norm(Perr)>1e-2)
        flag_err=1;
    end
end