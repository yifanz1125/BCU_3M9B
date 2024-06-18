%close all
%clear
%Cal_MM_Static;
%clear pfdata EMF Yload Case netdata

omegab=Basevalue.omegab;
%% Tfault and Tclear set
    Iter.Tfault=20;
    Iter.Trecover=20.247;
    Iter.Ttotal=100;
    Iter.Tunit=1e-4;
%% iteration procedure
% prefault
    delta0=prefault.SEP_delta;%zeros(ngen,1)D1D
    omega0=prefault.SEP_omegapu*omegab*ones(ngen,1);
    [theta_pre,omega_pre,thetac_pre,omegacoi_pre,Pe_pre,cycle_pre]=Fun_TrajIter_SRF(Iter.Tfault,Iter.Tunit,prefault.Yred,preset,delta0,omega0,omegab);
% faults
    delta0=theta_pre(cycle_pre,:);
    omega0=omega_pre(cycle_pre,:);
    [theta_fault,omega_fault,thetac_fault,omegacoi_fault,Pe_fault,cycle_fault]=Fun_TrajIter_SRF(Iter.Trecover-Iter.Tfault,Iter.Tunit,fault.Yred,preset,delta0,omega0,omegab);
% postfault
    delta0=theta_fault(cycle_fault,:);
    omega0=omega_fault(cycle_fault,:);
    [theta_post,omega_post,thetac_post,omegacoi_post,Pe_post,cycle_post]=Fun_TrajIter_SRF(Iter.Ttotal-Iter.Trecover,Iter.Tunit,postfault.Yred,preset,delta0,omega0,omegab);
    clear delta0 omega0

%% Data Collection
    Tpre=Iter.Tfault;
    Tfault=Iter.Trecover-Iter.Tfault;
    Tpost=Iter.Ttotal-Iter.Trecover;
    Tunit=Iter.Tunit;
    Ttotal=Iter.Ttotal;
    TM_pre=0:Tunit:Tpre-Tunit;
    TM_fault=Tpre:Tunit:Tpre+Tfault-Tunit;
    TM_post=round((Tpre+Tfault)/Tunit)*Tunit:Tunit:(round((Ttotal)/Tunit)-1)*Tunit;
    clear Tpre Tfault Tpost Tunit
    n_tt=cycle_pre+cycle_fault+cycle_post;
    IterData.omega=zeros(n_tt,ngen);
    IterData.omegacoi=zeros(n_tt,1);
    IterData.omegac=zeros(n_tt,ngen);
    IterData.theta=zeros(n_tt,ngen);
    IterData.thetac=zeros(n_tt,ngen);
    IterData.Tout=zeros(n_tt,1);
    IterData.Pe=zeros(n_tt,ngen);
    IterData.theta(1:cycle_pre,:)=theta_pre;
    IterData.thetac(1:cycle_pre,:)=thetac_pre;
    IterData.omega(1:cycle_pre,:)=omega_pre/omegab;
    IterData.omegacoi(1:cycle_pre)=omegacoi_pre/omegab;
    for i=1:ngen
        IterData.omegac(1:cycle_pre,i)=omega_pre(:,i)/omegab-omegacoi_pre/omegab;
    end
    IterData.Tout(1:cycle_pre)=TM_pre;
    IterData.Pe(1:cycle_pre,:)=Pe_pre;

    IterData.theta(cycle_pre+1:cycle_pre+cycle_fault,:)=theta_fault;
    IterData.thetac(cycle_pre+1:cycle_pre+cycle_fault,:)=thetac_fault;
    IterData.omega(cycle_pre+1:cycle_pre+cycle_fault,:)=omega_fault/omegab;
    IterData.omegacoi(cycle_pre+1:cycle_pre+cycle_fault)=omegacoi_fault/omegab;
    for i=1:ngen
        IterData.omegac(cycle_pre+1:cycle_pre+cycle_fault,i)=omega_fault(:,i)/omegab-omegacoi_fault/omegab;
    end
    IterData.Tout(cycle_pre+1:cycle_pre+cycle_fault)=TM_fault;
    IterData.Pe(cycle_pre+1:cycle_pre+cycle_fault,:)=Pe_fault;

    IterData.theta(cycle_pre+cycle_fault+1:n_tt,:)=theta_post;
    IterData.thetac(cycle_pre+cycle_fault+1:n_tt,:)=thetac_post;
    IterData.omega(cycle_pre+cycle_fault+1:n_tt,:)=omega_post/omegab;
    IterData.omegacoi(cycle_pre+cycle_fault+1:n_tt)=omegacoi_post/omegab;
    for i=1:ngen
        IterData.omegac(cycle_pre+cycle_fault+1:n_tt,i)=omega_post(:,i)/omegab-omegacoi_post/omegab;
    end
    IterData.omegacd=IterData.omegac*preset.d/sum(preset.d);
    IterData.deltacd=IterData.thetac*preset.d/sum(preset.d);
    
    IterData.Tout(cycle_pre+cycle_fault+1:n_tt)=TM_post;
    IterData.Pe(cycle_pre+cycle_fault+1:n_tt,:)=Pe_post;
    clear n_tt
%     % relative theta
%     IterData.theta12=IterData.theta(:,1)-IterData.theta(:,2);
%     IterData.theta13=IterData.theta(:,1)-IterData.theta(:,3);
%     IterData.theta23=IterData.theta(:,2)-IterData.theta(:,3);

    %% Trajectories in COI
    % thetac
        figure(3);
        set(gca,'position',[0.115,0.12,0.815,0.84]);
        set(gcf,'position',[60 200 600 450]);
%         for i=1:ngen
%             plot(IterData.Tout,IterData.thetac(:,i),'linewidth',2,'color',[(50+20*i)/255 150/255 (250-20*i)/255]);    hold on;
%     %         plot(IterData.Tout,IterData.thetac(:,2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%     %         plot(IterData.Tout,IterData.thetac(:,3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
% %             xlim([19,25]);
%             strlab(i)="M"+i;
%             legend(strlab);
%         end
        n_start=cycle_pre-fix(10/Iter.Tunit)+1;
        n_end=cycle_pre+fix(40/Iter.Tunit);
        
        plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.thetac((n_start:n_end),3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
        %plot(IterData.Tout(n_start:n_end),IterData.deltacd((n_start:n_end),1),'linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
        
        grid on;    grid minor;
        ylim([-2,4]);
        % xlabel ylabel position
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        ylab_x=xmin-(xmax-xmin)/15;
        ylab_y=(ymax+ymin)/2;
        xlab_x=(xmax+xmin)/2;
        xlab_y=ymin-(ymax-ymin)/15;
        % during-fault area identification
        trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
        % axis font
        ax=gca;
        ax.FontName='Arial';
        ax.FontSize=14;
        xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
        ylabel('\delta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
%         legend('\theta_1','\theta_2','\theta_3');

    % omega
%         figure(2);
%         set(gca,'position',[0.115,0.12,0.815,0.84]);
%         set(gcf,'position',[60 200 600 450]);
%         for i=1:ngen
%             plot(IterData.Tout,IterData.omega(:,i),'linewidth',2,'color',[(50+20*i)/255 150/255 (250-20*i)/255]);    hold on;
%     %         plot(IterData.Tout,IterData.omega(:,2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%     %         plot(IterData.Tout,IterData.omega(:,3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
%             strlab(i)="M"+i;
%             legend(strlab);
%         end
%         grid on;    grid minor;
%         % xlabel ylabel position
% %         xlim([19,25]);
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         xl=xlim;
%         xmin=xl(1,1);
%         xmax=xl(1,2);
%         ylab_x=xmin-(xmax-xmin)/12;
%         ylab_y=(ymax+ymin)/2;
%         xlab_x=(xmax+xmin)/2;
%         xlab_y=ymin-(ymax-ymin)/15;
%         % during-fault area identification
%         trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         ylim([ymin,ymax]);
%         % axis font
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize=14;
%         xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
%         ylabel('\omega(pu)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
% %         legend('\omega_1','\omega_2','\omega_3');

%     % theta
%         figure(3);
%         set(gca,'position',[0.115,0.12,0.815,0.84]);
%         set(gcf,'position',[60 200 600 450]);
%         plot(IterData.Tout,IterData.theta(:,1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
%         plot(IterData.Tout,IterData.theta(:,2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
%         plot(IterData.Tout,IterData.theta(:,3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
% %         xlim([19,23]);
%         grid on;    grid minor;
%         % xlabel ylabel position
%         yl=ylim;
%         ymin=yl(1,1);
%         ymax=yl(1,2);
%         xl=xlim;
%         xmin=xl(1,1);
%         xmax=xl(1,2);
%         ylab_x=xmin-(xmax-xmin)/15;
%         ylab_y=(ymax+ymin)/2;
%         xlab_x=(xmax+xmin)/2;
%         xlab_y=ymin-(ymax-ymin)/15;
%         % during-fault area identification
%         trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
%         fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
%         % axis font
%         ax=gca;
%         ax.FontName='Arial';
%         ax.FontSize=14;
%         xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
%         ylabel('\theta_c(rad)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
%         legend('\theta_1','\theta_2','\theta_3');

        figure(4);
        set(gca,'position',[0.115,0.12,0.815,0.84]);
        set(gcf,'position',[60 200 600 450]);
        plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
        plot(IterData.Tout(n_start:n_end),IterData.omegac((n_start:n_end),3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
        
        plot(IterData.Tout(n_start:n_end),IterData.omegacoi((n_start:n_end),1)-1,'linewidth',2,'color',[0/255 0/255 0/255]);    hold on;
%         plot(IterData.Tout(n_start:n_end),IterData.omegacd((n_start:n_end),1),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
        grid on;    grid minor;
        ylim([-0.03,0.05]);
        % xlabel ylabel position
%         xlim([19,25]);
        yl=ylim;
        ymin=yl(1,1);
        ymax=yl(1,2);
        xl=xlim;
        xmin=xl(1,1);
        xmax=xl(1,2);
        ylab_x=xmin-(xmax-xmin)/12;
        ylab_y=(ymax+ymin)/2;
        xlab_x=(xmax+xmin)/2;
        xlab_y=ymin-(ymax-ymin)/15;
        % during-fault area identification
        trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];   thetarange=[ymin,ymin,ymax,ymax];
        fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5); hold on;
        ylim([ymin,ymax]);
        % axis font
        ax=gca;
        ax.FontName='Arial';
        ax.FontSize=14;
        xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
        ylabel('\omega(pu)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);

        
