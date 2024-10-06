%% Function: update start point along the ray of deltac_last and postfault SEP
function [deltac_update,theta_update,voltage_update,flag_update]=Fun_Cal_UpdateStartPoint_SPM(deltac_lastpoint,theta_lastpoint,voltage_lastpoint,preset,postfault)
%% Settings
    Yfull=postfault.Yfull_mod;
    ngen=preset.ngen;
    nbus=preset.nbus;
    G_post=real(Yfull);
    B_post=imag(Yfull);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_update=0;
    Basevalue = evalin('base','Basevalue');
%% initialization
    Pe=zeros(ngen,1);
    len_ray=1e-3;
    dir_ray=zeros(1,ngen);  % vector from SEP_post to deltac_lastpoint
    dir_ray_theta=zeros(1,nbus); 
    dir_ray_voltage=zeros(1,nbus);  

    deltac_SEP=postfault.SEP_delta;
    theta_SEP=postfault.net_delta;
    voltage_SEP=postfault.net_voltage;
    dir_ray=deltac_lastpoint-deltac_SEP;
    dir_ray_theta=theta_lastpoint-theta_SEP;
    dir_ray_voltage=voltage_lastpoint-voltage_SEP;

    n_itermax=fix(2*norm(dir_ray)/len_ray);
    len_ray_theta = 2*norm(dir_ray_theta)/n_itermax;
    len_ray_voltage = 2*norm(dir_ray_voltage)/n_itermax;
    dir_ray=dir_ray/norm(dir_ray);
    dir_ray_theta=dir_ray_theta/norm(dir_ray_theta);
    dir_ray_voltage=dir_ray_voltage/norm(dir_ray_voltage);

    n_iter=1;
    Ep_obsv=[0 0 0];

    Ep0=zeros(n_itermax,1);
    flag_position=0;    %1--init point is inside of boundary 2--outside

    for  i=1:ngen
        for j=1:ngen
            ddelta=deltac_lastpoint(i)-deltac_lastpoint(j);
            Pe(i)=Pe(i)+E(i)*E(j)*(G_post(i,j)*cos(ddelta)+B_post(i,j)*sin(ddelta));
        end
        for l=1:(nbus-ngen)
            ddelta=deltac_lastpoint(i)-theta_lastpoint(l);
            Pe(i)=Pe(i)+E(i)*voltage_lastpoint(l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_lastpoint(l)*G_post(i,l+ngen)*cos(ddelta);
        end
    end
    
%% Check the direction (Use current point's vector)
    dotproduct=dir_ray'*(Pm-Pe-m/sum(m)*sum(Pm-Pe));
    if(dotproduct>0)
        flag_position=2;
%         len_ray=-len_ray;
    else
        flag_position=1;
    end
%% Search for the local maximum Ep along ray(and extension)
    [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_SEP,theta_SEP,voltage_SEP);
    Ep_obsv(2)=sum(Ep);
    deltac_act=deltac_SEP+len_ray*dir_ray;
    theta_est = theta_SEP;
    voltage_est = voltage_SEP;
    [theta_act,voltage_act,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_act,preset,Basevalue,"postfault",1e4,1e-10);

    [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_act,theta_act,voltage_act);
    Ep_obsv(3)=sum(Ep);
    
    while(n_iter~=-1)
        Ep_obsv(1)=Ep_obsv(2);
        Ep_obsv(2)=Ep_obsv(3);
        deltac_act=deltac_act+len_ray*dir_ray;
        theta_est = theta_act;
        voltage_est = voltage_act;
        [theta_act,voltage_act,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_act,preset,Basevalue,"postfault",1e4,1e-10);

        [Ep(1),Ep(2),Ep(3),Ep(4),Ep(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_act,theta_act,voltage_act);
        Ep_obsv(3)=sum(Ep);
        Ep0(n_iter)=Ep_obsv(3);
        n_iter=n_iter+1;
        if(Ep_obsv(2)>Ep_obsv(1)&&Ep_obsv(2)>Ep_obsv(3))
            n_iter=-1;
        elseif(n_iter>n_itermax)
            break;
        end
    end

    if(n_iter==-1)
        flag_update=1;
        deltac_update=deltac_act-len_ray*dir_ray;
        theta_est = theta_act;
        voltage_est = voltage_act;
        [theta_update,voltage_update,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_update,preset,Basevalue,"postfault",1e4,1e-10);
    else
%         error('No local maximum point found!');
        flag_update=0;
        deltac_update=deltac_lastpoint-len_ray*dir_ray;
        [theta_update,voltage_update,flag_iter,n_it,err] = Fun_AEiteration_SPM(theta_est,voltage_est,deltac_update,preset,Basevalue,"postfault",1e4,1e-10);
    end
end
