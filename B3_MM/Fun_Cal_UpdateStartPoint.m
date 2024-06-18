%% Function: update start point along the ray of thetac_last and postfault SEP
function [thetac_update,flag_update]=Fun_Cal_UpdateStartPoint(thetac_lastpoint,preset,postfault)
%% Settings
    Y_red=postfault.Yred;
    ngen=size(Y_red,1);
    G_post=real(Y_red);
    B_post=imag(Y_red);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_update=0;
%% initialization
    Pe=zeros(ngen,1);
    len_ray=1e-3;
    dir_ray=zeros(1,ngen);  % vector from SEP_post to thetac_lastpoint
    thetac_SEP=postfault.SEP_delta;
    dir_ray=thetac_lastpoint-thetac_SEP;
    n_itermax=fix(2*norm(dir_ray)/len_ray);
    dir_ray=dir_ray/norm(dir_ray);
    n_iter=1;
    Ep_obsv=[0 0 0];

    Ep0=zeros(n_itermax,1);
    flag_position=0;    %1--init point is inside of boundary 2--outside
    for i=1:ngen
        for j=1:ngen
            ddelta=thetac_lastpoint(i)-thetac_lastpoint(j);
            Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
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
    [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_SEP,thetac_SEP);
    Ep_obsv(2)=sum(Ep);
    thetac_act=thetac_SEP+len_ray*dir_ray;
%     thetac_act=thetac_SEP+len_ray*dir_ray;
    [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_SEP,thetac_act);
    Ep_obsv(3)=sum(Ep);
    
    while(n_iter~=-1)
        Ep_obsv(1)=Ep_obsv(2);
        Ep_obsv(2)=Ep_obsv(3);
        thetac_act=thetac_act+len_ray*dir_ray;
        [Ep(1),Ep(2),Ep(3)]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_SEP,thetac_act);
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
        thetac_update=thetac_act-len_ray*dir_ray;
    else
%         error('No local maximum point found!');
        flag_update=0;
        thetac_update=thetac_lastpoint-len_ray*dir_ray;
    end
end
