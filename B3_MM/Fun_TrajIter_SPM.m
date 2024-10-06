function [delta_RK4,omega_RK4,deltac_RK4,omegacoi_RK4,theta_net_RK4,voltage_net_RK4,cycle]=Fun_TrajIter_SPM(Tlength,Tunit,system,preset,delta0,omega0,theta_net0,voltage_net0,Basevalue)
    cycle=round(Tlength/Tunit);
    ngen=preset.ngen;
    nbus=preset.nbus;
    omega_RK4=zeros(cycle,ngen);
    theta_net_RK4=zeros(cycle,nbus-ngen);
    voltage_net_RK4=zeros(cycle,nbus-ngen);
    omegac_RK4=zeros(cycle,ngen);
    delta_RK4=zeros(cycle,ngen);
    deltac_RK4=zeros(cycle,ngen);
    omegacoi_RK4=zeros(cycle,1);
    m=preset.m;
    mT=sum(m,1);
    K_deltaomega=zeros(ngen*2,4);
    RK4ratio=[1;2;2;1]/6;
    flag_wrap=zeros(ngen,1);
%% initialize
    delta_RK4(1,:)=delta0;
    omega_RK4(1,:)=omega0;
    theta_net_RK4(1,:)=theta_net0;
    voltage_net_RK4(1,:)=voltage_net0;

%% Iteration procedure
    for tm=1:cycle
    % current state
    net_value = [theta_net_RK4(tm,:), voltage_net_RK4(tm,:)]';
    x = [delta_RK4(tm,:), omega_RK4(tm,:)]';

    % Runge-Kutta 4th order
    K_deltaomega(:,1)=F_3M9B_SP_ODE(x,net_value,system)*Tunit;
    K_deltaomega(:,2)=F_3M9B_SP_ODE(x+K_deltaomega(:,1)/2,net_value,system)*Tunit;
    K_deltaomega(:,3)=F_3M9B_SP_ODE(x+K_deltaomega(:,2)/2,net_value,system)*Tunit;
    K_deltaomega(:,4)=F_3M9B_SP_ODE(x+K_deltaomega(:,3),net_value,system)*Tunit;
   
    %update x
    x_new = x + K_deltaomega*RK4ratio;
    x_new=x_new';

    if(tm<cycle)
       delta_RK4(tm+1,:) = x_new(1:3);
       omega_RK4(tm+1,:) = x_new(4:6);
       % update net_value
       [net_value_new,fval,exitflag,output]=fsolve(@(x)Fun_AEfslove_SPM(x,delta_RK4(tm+1,:)',preset,system),net_value,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
       if exitflag<=0
               disp('fval=');
               disp(fval);
               disp('tm=');
               disp(tm);
               disp('error=');
               disp(maxabs(fval));
               error('cannot find fault-clear state! \n');
       else
           theta_net_RK4(tm+1,:) = net_value_new(1:6)';
           voltage_net_RK4(tm+1,:) = net_value_new(7:12)';
       end
              
%        [delta_net_s,voltage_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(net_value(1:6),net_value(7:12),delta_RK4(tm+1,:),preset,Basevalue,system,1e4,1e-10);
%        net_value_new = [delta_net_s;voltage_net_s];
%        if flag_iter~=1
%                disp('flag_iter=');
%                disp(flag_iter);
%                disp('tm=');
%                disp(tm);
%                disp('error=');
%                disp(maxabs(err));
%                error('cannot find fault-clear state! \n');
%        else
%            theta_net_RK4(tm+1,:) = net_value_new(1:6)';
%            voltage_net_RK4(tm+1,:) = net_value_new(7:12)';
%        end
    end
    
    % theta wrap procedure
        for i=1:ngen
            if(delta_RK4(tm,i)>pi&&flag_wrap(i)==0)    
                flag_wrap(i)=1;
            end
        end
        if(isequal(flag_wrap,ones(ngen,1))==1)
            delta_RK4(tm,:)=delta_RK4(tm,:)-2*pi*ones(1,ngen);
            if(tm<cycle)
                delta_RK4(tm+1,:)=delta_RK4(tm+1,:)-2*pi*ones(1,ngen);
            end
            flag_wrap=zeros(ngen,1);
        end   
        omegacoi_RK4(tm,1)=omega_RK4(tm,:)*m/mT;
        deltacoi=delta_RK4(tm,:)*m/mT;
        deltac_RK4(tm,:)=delta_RK4(tm,:)-deltacoi*ones(1,ngen);
        omegac_RK4(tm,:)=omega_RK4(tm,:)-omegacoi_RK4(tm,1)*ones(1,ngen);

    end
end
