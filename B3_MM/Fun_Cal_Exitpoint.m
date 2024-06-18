function [theta_RK4,omega_RK4,thetac_RK4,omegac_RK4,Exittm]=Fun_Cal_Exitpoint(Tlength,Tunit,Yredfault,Yredpost,preset,delta0,omega0,omegab)
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_exit=0;
    cycle=fix(Tlength/Tunit);
    G_RK4=real(Yredfault);
    B_RK4=imag(Yredfault);
    G_post=real(Yredpost);
    B_post=imag(Yredpost);
    ngen=size(Yredfault,1);
    Pe=zeros(cycle,ngen);
    Pe_post=zeros(cycle,ngen);
    omega_RK4=zeros(cycle,ngen);
    omegac_RK4=zeros(cycle,ngen);
    theta_RK4=zeros(cycle,ngen);
    thetac_RK4=zeros(cycle,ngen);
    omegacoi_RK4=zeros(cycle,1);
    mT=sum(m,1);
    S_=zeros(ngen,4);
    S=zeros(ngen,4);
    RK4ratio=[1;2;2;1]*Tunit/6;
    flag_wrap=zeros(ngen,1);
    Dotproduct=[0 0];
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
        omegacoi_RK4(tm,1)=omega_RK4(tm,:)*m/mT;
        thetacoi=theta_RK4(tm,:)*m/mT;
        thetac_RK4(tm,:)=theta_RK4(tm,:)-thetacoi*ones(1,ngen);
        omegac_RK4(tm,:)=omega_RK4(tm,:)-omegacoi_RK4(tm,1)*ones(1,ngen);
    % theta wrap procedure
        for i=1:ngen
            if(theta_RK4(tm,i)>pi&&flag_wrap(i)==0)    
                flag_wrap(i)=1;
            end
        end
        if(isequal(flag_wrap,ones(ngen,1))==1)
            theta_RK4(tm,:)=theta_RK4(tm,:)-2*pi*ones(1,ngen);
            if(tm<cycle)
                theta_RK4(tm+1,:)=theta_RK4(tm+1,:)-2*pi*ones(1,ngen);
            end
            flag_wrap=zeros(ngen,1);
        end   
    

        % Pe_post calculation
            for i=1:ngen
                for j=1:ngen
                    delta=theta_RK4(tm,i)-theta_RK4(tm,j);
                    Pe_post(tm,i)=Pe_post(tm,i)+E(i)*E(j)*(G_post(i,j)*cos(delta)+B_post(i,j)*sin(delta));
                end
            end
        Dotproduct(1)=Dotproduct(2);
        Dotproduct(2)=(Pm'-Pe_post(tm,:))*omegac_RK4(tm,:)';
        if(flag_exit==0&&tm>2&&Dotproduct(1)<0&&Dotproduct(2)>0)
            flag_exit=1;
            Exittm=tm-1;
        end
        
%             Dotproduct(tm+1)=(Pm'-Pe_post(tm+1,:))*omegac_RK4(tm+1,:)';
    end
end
