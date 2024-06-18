%% function: Calculate Damping energy and (Ek,Ep) of critical trajectory
%% preorder code: Fun_Cal_CCT_Real
%% output: Ed1--uniform damping energy; Ed2--nonuniform damping energy
function [Energy,Group]=Fun_Cal_DampingEnergy(postfault,preset,Basevalue,Critical)
%% Configurations
    ngen=size(postfault.Yred,1);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    Yred_post=postfault.Yred;
    G_post=real(Yred_post);
    B_post=imag(Yred_post);
    omegab=Basevalue.omegab;
    if(mod(ngen,2)==1)
        flag_odd=1; % number of generators is odd
    else
        flag_odd=0;
    end
    cycle_stb=size(Critical.Traj.Stb.thetac,1);
    cycle_unstb=size(Critical.Traj.Unstb.thetac,1);
%% active length of trajectories (Tunit=1e-4)
    cycle_stb=5e4;
    cycle_unstb=5e4;

    Critical.Traj.Stb.thetac=Critical.Traj.Stb.thetac(1:cycle_stb,:);
    Critical.Traj.Stb.omegac=Critical.Traj.Stb.omegac(1:cycle_stb,:);
    Critical.Traj.Stb.theta=Critical.Traj.Stb.theta(1:cycle_stb,:);
    Critical.Traj.Stb.omega=Critical.Traj.Stb.omega(1:cycle_stb,:);
    Critical.Traj.Stb.omegacoi=Critical.Traj.Stb.omega(1:cycle_stb,:)*m/sum(m);
    Critical.Traj.Unstb.thetac=Critical.Traj.Unstb.thetac(1:cycle_unstb,:);
    Critical.Traj.Unstb.omegac=Critical.Traj.Unstb.omegac(1:cycle_unstb,:);
    Critical.Traj.Unstb.theta=Critical.Traj.Unstb.theta(1:cycle_unstb,:);
    Critical.Traj.Unstb.omega=Critical.Traj.Unstb.omega(1:cycle_unstb,:);
%% Find a point on critical-traj closed to CUEP
    %% for ngen, candidate groups numbers=2^(ngen-1)-1
    n_group=2^(ngen-1)-1;
    Group.delta_stb=zeros(cycle_stb,n_group);
    Group.delta_unstb=zeros(cycle_unstb,n_group);
    Group.omega_stb=zeros(cycle_stb,n_group);
    Group.omega_unstb=zeros(cycle_unstb,n_group);
    
    Group.CUEP=zeros(1,n_group);
    Group.gennum=zeros(1,n_group);
    Group.genno=size(fix(ngen/2),n_group);
    Group.flag_reverse=zeros(1,n_group);    % 0--[Group.genno] stands for S-group, 1--[Group.genno] stands for A-group

    %% All candidate groups
    no_gen=1:1:ngen;
    no_group=1;    % no of group (cosistent with delta_group)
    if(flag_odd==1)
        num_layer=(ngen-1)/2;
    else
        num_layer=ngen/2;
    end
    % overwrite layer: select specific groups
%     num_layer=2;

        for n_layer=1:num_layer
            Currentlayer.n_grp=nchoosek(ngen,n_layer);
            Currentlayer.candidategrp=nchoosek(no_gen,n_layer);
            if(flag_odd==0&&n_layer==num_layer) % deal with C(n/2)(n)
                Currentlayer.candidategrp(Currentlayer.n_grp/2+1:Currentlayer.n_grp,:)=[];
                Currentlayer.n_grp=Currentlayer.n_grp/2;
            end
            for no_ingrp=1:Currentlayer.n_grp
                Currentlayer.currentgrp=Currentlayer.candidategrp(no_ingrp,:);
                m_g2=preset.m;
                for n_slt=1:size(Currentlayer.currentgrp,2)
                    m_g2(Currentlayer.currentgrp(1,n_slt))=0;
                    Group.genno(n_slt,no_group)=Currentlayer.currentgrp(1,n_slt);
                end
                m_g1=preset.m-m_g2;
                % calculate current group's delta
                Group.delta_stb(:,no_group)=Critical.Traj.Stb.thetac*m_g1/sum(m_g1)-Critical.Traj.Stb.thetac*m_g2/sum(m_g2);
                Group.delta_unstb(:,no_group)=Critical.Traj.Unstb.thetac*m_g1/sum(m_g1)-Critical.Traj.Unstb.thetac*m_g2/sum(m_g2);
                
                Group.omega_stb(:,no_group)=Critical.Traj.Stb.omegac*m_g1/sum(m_g1)-Critical.Traj.Stb.omegac*m_g2/sum(m_g2);
                Group.omega_unstb(:,no_group)=Critical.Traj.Unstb.omegac*m_g1/sum(m_g1)-Critical.Traj.Unstb.omegac*m_g2/sum(m_g2);
                
                Group.CUEP(no_group)=postfault.CUEP_delta'*m_g1/sum(m_g1)-postfault.CUEP_delta'*m_g2/sum(m_g2);
                if(Group.CUEP(no_group)<0)
                    Group.CUEP(no_group)=-Group.CUEP(no_group);
                    Group.delta_stb(:,no_group)=-Group.delta_stb(:,no_group);
                    Group.delta_unstb(:,no_group)=-Group.delta_unstb(:,no_group);
                    
                    Group.omega_stb(:,no_group)=-Group.omega_stb(:,no_group);
                    Group.omega_unstb(:,no_group)=-Group.omega_unstb(:,no_group);
                    Group.flag_reverse(no_group)=1;
                end
                Group.gennum(no_group)=n_layer;
                no_group=no_group+1;
            end
        end
        clear flag_odd num_layer n_layer no_ingrp n_slt m_g1 m_g2
        cnt_act=1;
%         no_CUEPact=zeros(size(CUEP,2),1);
        no_CUEPact=[];
        for i=1:size(Group.CUEP,2)
            if(Group.CUEP(i)>pi/2)
                no_CUEPact(cnt_act)=i;
                cnt_act=cnt_act+1;
            end
        end
        if(isempty(no_CUEPact)==0)
            Group.CUEPact=no_CUEPact;
        else
            Group.CUEPact=[];
        end
        clear no_CUEPact
        delta_stbmax=zeros(n_group,1);
        for i=1:n_group
            delta_stbmax(i)=max(Group.delta_stb(:,i));
        end
        Group.deltamax_stb=delta_stbmax;
        clear delta_stbmax

    
    %% select critical group
    % for 3M
%         Group.omegainitlast_unstb=zeros(1,n_group);
%         Group.omegainitlast_unstb(1)=abs(Critical.Traj.Unstb.omegac(cycle_unstb,2)-Critical.Traj.Unstb.omegac(cycle_unstb,3));
%         Group.omegainitlast_unstb(2)=abs(Critical.Traj.Unstb.omegac(cycle_unstb,1)-Critical.Traj.Unstb.omegac(cycle_unstb,3));
%         Group.omegainitlast_unstb(3)=abs(Critical.Traj.Unstb.omegac(cycle_unstb,1)-Critical.Traj.Unstb.omegac(cycle_unstb,2));
%         [num,crt_grp]=min(Group.omegainitlast_unstb);
%         Group.criticalgroup=crt_grp;
%         clear num
    % for 10M
    for no_group=1:3
        figure(1);
        plot(Group.delta_stb(:,no_group)); hold on;
        plot(Group.CUEP(no_group)*ones(cycle_stb,1),':');    hold on;
        strlb(2*no_group-1)="group"+no_group;
        strlb(2*no_group)="CUEP"+no_group;
        legend(strlb);        
        figure(3)
        plot(Group.delta_unstb(:,no_group)); hold on;
        plot(Group.CUEP(no_group)*ones(cycle_stb,1),':');    hold on;
        strlb(no_group*2-1)="group"+no_group;
        strlb(no_group*2)="CUEP"+no_group;
        legend(strlb);
    end
    if(isempty(Group.CUEPact)==0)
%     for i=1:size(Group.CUEPact,2)
%         figure(1);
%         plot(Group.delta_stb(:,Group.CUEPact(i)),'linewidth',2); hold on;
%         plot(Group.CUEP(Group.CUEPact(i))*ones(cycle_stb,1),'linewidth',2); hold on;
%         figure(3);
%         plot(Group.delta_unstb(:,Group.CUEPact(i)),'linewidth',2); hold on;
%         plot(Group.CUEP(Group.CUEPact(i))*ones(cycle_stb,1),'linewidth',2); hold on;
%     end
       fprintf('Groups in which UEP>pi/2 are: ');
       fprintf('%d ',Group.CUEPact);
    end
       crt_grp=0;
       fprintf('\n');
       while(crt_grp<1||crt_grp>n_group)
           crt_grp=input('Critical group is selected as:');
           if(crt_grp<1||crt_grp>n_group)
               fprintf('Input data is invalid!\n');
               crt_grp=0;
           end
       end

    Group.criticalgroup=crt_grp;
    %% find end point along critical group
        dis_min=Group.CUEP(crt_grp);
        no_omegazero=0;
        no_overCUEP=-1;
        no_max=0; % maximum delta @first swing
        no_closeUEP=0;  % closed to CUEP (positive--delta doesn't exceed CUEP, negative--delta exceeds CUEP)
        dis=Group.CUEP(crt_grp)-Group.delta_stb(1,crt_grp);
        norm_fulldimension=zeros(cycle_stb,1);
        norm_fulldimension(1)=norm(postfault.CUEP_delta-Critical.Traj.Stb.thetac(1,:));
        for tm=2:cycle_stb
            % search for omega_group cross zero point
            if(Group.omega_stb(tm-1,crt_grp)>0&&Group.omega_stb(tm,crt_grp)<0&&no_omegazero==0)
                no_omegazero=tm-1;
                no_closeUEP=tm-1;       % for traj doesn't exceed CUEP
            end
            dis_his=dis;
            dis=Group.CUEP(crt_grp)-Group.delta_stb(tm,crt_grp);
            norm_fulldimension(tm)=norm(postfault.CUEP_delta-Critical.Traj.Stb.thetac(tm,:));
            if(dis<dis_min)
                dis_min=dis;
                if(dis_min<0&&no_overCUEP==-1)
                    no_overCUEP=tm;
                end
            end
            if(dis>dis_his&&no_max==0)
                no_max=tm-1;
            end
            if(no_omegazero~=-1&&abs(dis)<abs(dis_his))
                no_closeUEP=tm; % for traj exceeds CUEP and turns back to CUEP
            end
            if(no_omegazero~=0&&no_max~=0&&no_closeUEP~=0&&dis>0)
                break;
            end
        end

        if(no_overCUEP==-1) % tractory didn't exceed CUEP
            no_critical=no_omegazero;
            Group.flag_overCUEP=0;
        else
%             no_critical=no_overCUEP;
            no_critical=no_omegazero;
%             no_critical=no_closeUEP;
            Group.flag_overCUEP=1;
        end


%         
%         clear dis dis_his crt_grp
        Critical.Traj.Stb.thetad=Critical.Traj.Stb.thetac*preset.d/sum(preset.d);   % nonuniform damping theta
        Critical.Traj.Stb.omegacoi=Critical.Traj.Stb.omega*preset.m/sum(preset.m);   % nonuniform damping theta
%         no_critical=1.7e4;
       fprintf('no_critical is ');
       fprintf('%d ',no_critical);
       fprintf('\n');
        no_overwrite=input('Overwrite critical point:');
        if(isempty(no_overwrite)==0)
            if(no_overwrite>0&&no_overwrite<cycle_stb)
                no_critical=no_overwrite;
            end
        end
%         clear no_overwrite
%         no_critical=1500;
%% Calculate damping energy along postfault traj
    %% Potential energy
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,Critical.Traj.Stb.thetac(1,:)');
        Energy.Ep_start=sum(Ep_tmp);
        Ep_lossy_start=Ep_tmp(3);
        Ep_mag_start=Ep_tmp(2);
        clear Ep_tmp
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,Critical.Traj.Stb.thetac(no_critical,:)');
        Energy.Ep_end=sum(Ep_tmp);   clear Ep_tmp
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,postfault.CUEP_delta);
        Energy.Ep_CUEP=sum(Ep_tmp);   clear Ep_tmp
        Energy.Ep=zeros(2e4,1);
        Energy.Ep_path=zeros(2e4,1);
        Energy.Ep_mag=zeros(2e4,1);
        for tm=1:size(Energy.Ep,1)
            [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3)]=Fun_Cal_PotentialEnergy(preset,postfault,postfault.SEP_delta,Critical.Traj.Stb.thetac(tm,:)');
            Energy.Ep_path(tm)=Ep_tmp(3);
            Energy.Ep_mag(tm)=Ep_tmp(2);
            Energy.Ep(tm)=sum(Ep_tmp);   clear Ep_tmp
        end
    %% Kinetic energy
        Ek0=zeros(ngen,1);
        Eke=zeros(ngen,1);
        for i=1:ngen
            Ek0(i)=0.5*m(i)*Critical.Traj.Stb.omegac(1,i)^2;
            Eke(i)=0.5*m(i)*Critical.Traj.Stb.omegac(no_critical,i)^2;
        end
        Energy.Ek_start=sum(Ek0);
        Energy.Ek_end=sum(Eke);
   %% Damping energy in iteration
        Ed_un=zeros(ngen,1);
        Ed_non=0;
        Ep_iter=zeros(ngen,1);
        % for Potential energy error display (in Fig5-7)
        Ep_lossy_iter=zeros(no_critical,1);
        Ep_mag_iter=zeros(no_critical,1);
        Ep_iter_record=zeros(no_critical,1);
        Ep_iter_record(1)=Energy.Ep_start;   
        Ep_lossy_iter(1)=Ep_lossy_start;
        Ep_mag_iter(1)=Ep_mag_start;
        Ek=zeros(no_critical,1);
        Ek(1)=Energy.Ek_start;
        Err_step=zeros(no_critical,1);

        for tm=2:no_critical
            ddelta_tmp=Critical.Traj.Stb.thetac(tm,:)-Critical.Traj.Stb.thetac(tm-1,:);
            Pe=zeros(ngen,1);
            Pe_lossy=zeros(ngen,1);  % lossy dispattive energy
            Pe_mag=zeros(ngen,1);   % Pe relevant to line inductance
            for i=1:ngen
                for j=1:ngen
                    delta=Critical.Traj.Stb.thetac(tm,i)-Critical.Traj.Stb.thetac(tm,j);
                    Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(delta)+E(i)*E(j)*G_post(i,j)*cos(delta);
                    if(i~=j)
                    Pe_lossy(i)=Pe_lossy(i)+E(i)*E(j)*G_post(i,j)*cos(delta);
                    Pe_mag(i)=Pe_mag(i)+E(i)*E(j)*B_post(i,j)*sin(delta);
                    end
                end
            end

            for i=1:ngen
                Ed_un(i)=Ed_un(i)+preset.d(i)*Critical.Traj.Stb.omegac(tm-1,i)*ddelta_tmp(i);
                Ed_non=Ed_non+preset.d(i)*(Critical.Traj.Stb.omegacoi(tm-1,1)-omegab)*ddelta_tmp(i);
                Ep_iter(i)=Ep_iter(i)-(Pm(i)-Pe(i))*ddelta_tmp(i);
                Ek(tm)=Ek(tm)+0.5*m(i)*Critical.Traj.Stb.omegac(tm,i)^2;
            end
            Ep_lossy_iter(tm)=Ep_lossy_iter(tm-1)+ddelta_tmp*Pe_lossy;
            Ep_mag_iter(tm)=Ep_mag_iter(tm-1)+ddelta_tmp*Pe_mag;          
            Ep_iter_record(tm)=Ep_iter_record(tm-1)-ddelta_tmp*(Pm-Pe);
          
            Err_step(tm)=sum(Ep_iter)+sum(Ed_un)+Ed_non+Ek(tm)-Ek(1);
        end



        Energy.Ep_iter=sum(Ep_iter);
        Energy.Ed_uniform=sum(Ed_un);
        Energy.Ed_nonuniform=Ed_non;


    
    figure(5);
    plot(Energy.Ep_path(1:no_critical),'color','r','LineWidth',2);  hold on;
    plot(Ep_lossy_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    legend('Potential','Iter');
    title('lossy energy');
    figure(6);
    plot(Energy.Ep_mag(1:no_critical),'color','r','LineWidth',2);  hold on;
    plot(Ep_mag_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    legend('Potential','Iter');
    title('magnetic energy');
    figure(7);
    plot(Energy.Ep_path(1:no_critical)-Ep_lossy_iter,'color','r','LineWidth',2);  hold on;
    plot(Energy.Ep(1:no_critical)-Ep_iter_record,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    legend('Lossy part err','Total err');
    title('Error');    

end


%% plot
%         figure(1);
%         plot(Group.delta_stb(:,Group.criticalgroup)); hold on;
%         plot(Group.CUEP(Group.criticalgroup)*ones(cycle_stb,1),':');    hold on;
%         strlb(Group.criticalgroup)="group"+Group.criticalgroup;
%         legend(strlb);
%         figure(2);
%         plot(Group.omega_stb(:,Group.criticalgroup)); hold on;
%         grid on;

%     for no_group=1:n_group
%         figure(1);
%         plot(Group.delta_stb(:,no_group)); hold on;
%         plot(Group.CUEP(no_group)*ones(cycle_stb,1),':');    hold on;
%         strlb(no_group)="group"+no_group;
%         legend(strlb);
%         figure(2);
%         plot(Group.omega_stb(:,no_group)); hold on;
%     end
    
        
    