function [theta_RK4,omega_RK4,thetac_RK4,omegacoi_RK4,Pe,cycle,flag_unstb]=Fun_TrajIter_StableCheck_SRF(Tlength,Tunit,postfault,preset,delta0,omega0,omegab)
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    cycle=round(Tlength/Tunit);
    Yred=postfault.Yred;
    SEP_post=postfault.SEP_delta;
    G_RK4=real(Yred);
    B_RK4=imag(Yred);
    ngen=size(Yred,1);
    Pe=zeros(cycle,ngen);
    omega_RK4=zeros(cycle,ngen);
    omegac_RK4=zeros(cycle,ngen);
    theta_RK4=zeros(cycle,ngen);
    thetac_RK4=zeros(cycle,ngen);
    omegacoi_RK4=zeros(cycle,1);
    mT=sum(m,1);
    S_=zeros(ngen,4);
    S=zeros(ngen,4);
    RK4ratio=[1;2;2;1]*Tunit/6;
%     flag_wrap=zeros(ngen,1);
    flag_unstb=0;
%% initialize
    theta_RK4(1,:)=delta0;
    omega_RK4(1,:)=omega0;
%% Iteration procedure
    for tm=1:cycle
    % Pe calculation
        for i=1:ngen
            for j=1:ngen
                delta=theta_RK4(tm,i)-theta_RK4(tm,j);
                Pe(tm,i)=Pe(tm,i)+E(i)*E(j)*(G_RK4(i,j)*cos(delta)+B_RK4(i,j)*sin(delta));
            end
        end
    % omega/theta derivation calculation
        for i=1:ngen
        % omega
            S_(i,1)=(Pm(i)-Pe(tm,i)-d(i)*(omega_RK4(tm,i)-omegab))/m(i);
            S_(i,2)=(Pm(i)-Pe(tm,i)-d(i)*(omega_RK4(tm,i)-omegab+S_(i,1)*Tunit/2))/m(i);
            S_(i,3)=(Pm(i)-Pe(tm,i)-d(i)*(omega_RK4(tm,i)-omegab+S_(i,2)*Tunit/2))/m(i);
            S_(i,4)=(Pm(i)-Pe(tm,i)-d(i)*(omega_RK4(tm,i)-omegab+S_(i,3)*Tunit))/m(i);
        % delta
            S(i,1)=omega_RK4(tm,i);%*omegab;
            S(i,2)=(omega_RK4(tm,i)+S_(i,1)*Tunit/2);%*omegab;
            S(i,3)=(omega_RK4(tm,i)+S_(i,2)*Tunit/2);%*omegab;
            S(i,4)=(omega_RK4(tm,i)+S_(i,3)*Tunit);%*omegab;
        end
    % omega/theta update
        if(tm<cycle)
            for i=1:ngen
                omega_RK4(tm+1,i)=omega_RK4(tm,i)+S_(i,:)*RK4ratio;
                theta_RK4(tm+1,i)=theta_RK4(tm,i)+S(i,:)*RK4ratio;
            end
        end
    % coi calculate
        omegacoi_RK4(tm,1)=omega_RK4(tm,:)*m/mT;
        thetacoi=theta_RK4(tm,:)*m/mT;
        thetac_RK4(tm,:)=theta_RK4(tm,:)-thetacoi*ones(1,ngen);
        omegac_RK4(tm,:)=omega_RK4(tm,:)-omegacoi_RK4(tm,1)*ones(1,ngen);
    % stability check
        deltamax=0;
        for i=1:ngen
            for j=1:ngen
                deltacurrent=theta_RK4(tm,i)-theta_RK4(tm,j);
                if(abs(deltacurrent)>deltamax)
                    deltamax=abs(deltacurrent);
                end
            end
        end
        if(deltamax>=2*pi)
            flag_unstb=1;
%             fprintf('Trajectory is judged unstable : Delta>2pi');
%             cycle=tm;
%             break;
        end
        % whether trajectory is closed to postfault_SEP
        if(tm==cycle&&flag_unstb==0)
            if(norm(thetac_RK4(tm,:)'-SEP_post)>0.1)  
                flag_unstb=1;
                fprintf('Trajectory is judged unstable: trajectory is apart from Postfault SEP @Tlength');
            end
        end

%     % theta wrap procedure
%         for i=1:ngen
%             if(theta_RK4(tm,i)>pi&&flag_wrap(i)==0)    
%                 flag_wrap(i)=1;
%             end
%         end
%         if(isequal(flag_wrap,ones(ngen,1))==1)
%             theta_RK4(tm,:)=theta_RK4(tm,:)-2*pi*ones(1,ngen);
%             if(tm<cycle)
%                 theta_RK4(tm+1,:)=theta_RK4(tm+1,:)-2*pi*ones(1,ngen);
%             end
%             flag_wrap=zeros(ngen,1);
%         end


    end
end
