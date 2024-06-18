%% function: calculate CCT along fault-on trajectory where Ek+Ep=Ecritical
function [CCT,Exit_thetac,Exit_omegac,Exit_theta,Exit_omega,flag_CCT]=Fun_Cal_CCT_Energy(E_critical,fault,postfault,preset)
    cycle=size(fault.traj.omega,1);
    Ek=zeros(cycle,1);
    Ep=zeros(cycle,1);
    Esum=zeros(cycle,1);
    ngen=size(postfault.Yred,1);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_CCT=0;
    for tm=1:cycle
        for i=1:ngen
        Ek(tm)=Ek(tm)+0.5*m(i)*fault.traj.omegac(tm,i)^2;
        end
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,fault.traj.thetac(tm,:)');
        Ep(tm)=sum(Ep_tmp);
        clear Ep_tmp
        Esum(tm)=Ep(tm)+Ek(tm);
        if(tm>1)
        if(Esum(tm-1)<E_critical&&Esum(tm)>E_critical&&flag_CCT==0)
            CCT=(tm-1)*fault.traj.Tunit;
            Exit_thetac=fault.traj.thetac(tm-1,:);
            Exit_omegac=fault.traj.omegac(tm-1,:);
            Exit_theta=fault.traj.theta(tm-1,:);
            Exit_omega=fault.traj.omega(tm-1,:);
            flag_CCT=1;
        end
        end
    end
end