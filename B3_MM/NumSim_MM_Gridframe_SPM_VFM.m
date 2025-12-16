%close all
%clear
%Cal_MM_Static;
%clear pfdata EMF Yload Case netdata

omegab=Basevalue.omegab;
%% Tfault and Tclear set
    Iter.Tfault=20;
    Iter.Trecover=20.21;%20.27;
    Iter.Ttotal=24;
    Iter.Tunit=1e-4;
    T_before = 5;
    T_after = 4; %10
    cycle_pre=round(Iter.Tfault/Iter.Tunit);
    n_start=cycle_pre-fix(T_before/Iter.Tunit)+1;
    n_end=cycle_pre+fix(T_after/Iter.Tunit);

%% Results derived by DAE function in Structure Preserved Model
M = diag([ones(3,1); ones(12,1)*1e-15; ones(3,1)]);
options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,3),1e-12*ones(1,12),1e-8*ones(1,3)]);
% prefault
    system="prefault";
    delta0=prefault.SEP_delta;
    omega0=prefault.SEP_omegapu*omegab*ones(2,1);
    delta_net0=prefault.net_delta;
    voltage_net0=prefault.net_voltage;
    x0=[delta0; delta_net0; voltage_net0; omega0(1);prefault.SEP_ydc;omega0(2)];
    [t_prefault, x_prefault_all] = ode15s(@(t,x)f_timedomain_DAE(t, x, system),0:Iter.Tunit:Iter.Tfault,x0,options);
    out2 = [];
    for i =2: length(t_prefault)
        out = F_3M9B_SP_DAE_VFM(x_prefault_all(i,:),system);
        out2 = [out2;out(2)];
    end
    omega_VFM_prefault = out2;

%% fault
    system="fault1";
    delta0=x_prefault_all(end,1:3)';
    ydc0 = x_prefault_all(end,17);
    omega0=[x_prefault_all(end,16)';x_prefault_all(end,18)']; %1,3
    delta_net0=x_prefault_all(end,4:9)';
    voltage_net0=x_prefault_all(end,10:15)';
    delta_net0(fault.faultbus-ngen)=[];
    voltage_net0(fault.faultbus-ngen)=[];
    delta_net0(6)=0;
    voltage_net0(6)=0;
    [delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(delta_net0,voltage_net0,delta0,preset,Basevalue,system,1e5,1e-10);
    [t_fault, x_fault_all] = ode15s(@(t, x) f_timedomain_DAE(t, x, system),Iter.Tfault:Iter.Tunit:Iter.Trecover,[delta0; delta_net_s; V_net_s; omega0(1);ydc0;omega0(2)],options);
    fault_delta_net=x_fault_all(2:end,4:8);
    fault_delta_net=[fault_delta_net(:,1:(fault.faultbus-ngen-1)) zeros(size(fault_delta_net,1),1) fault_delta_net(:,(fault.faultbus-ngen):end)];
    fault_voltage_net=x_fault_all(2:end,10:14);
    fault_voltage_net=[fault_voltage_net(:,1:(fault.faultbus-ngen-1)) zeros(size(fault_voltage_net,1),1) fault_voltage_net(:,(fault.faultbus-ngen):end)];
    fault_delta_gen = x_fault_all(2:end,1:3);
    out2 = [];
    for i =2: length(t_fault)
        out = F_3M9B_SP_DAE_VFM(x_fault_all(i,:),system);
        out2 = [out2;out(2)];
    end
    omega_VFM_fault = out2;

%% fault-clear
    system = "postfault";
    delta_net_faultclear = zeros(size(fault_delta_gen,1),nbus-ngen);
    voltage_net_faultclear = zeros(size(fault_delta_gen,1),nbus-ngen);
    delta_net_faultclear2 = zeros(size(fault_delta_gen,1),nbus-ngen);
    voltage_net_faultclear2 = zeros(size(fault_delta_gen,1),nbus-ngen);
    temp_ini = [zeros(6,1);ones(6,1)];
    for i = 1:(size(fault_delta_gen,1))
       %[delta_net_temp,voltage_temp,flag_iter,n_iter,err] = Fun_AEiteration_SPM(zeros(6,1),ones(6,1),fault_delta_gen(i,:)',preset,Basevalue,system,1e4,1e-10);
       %delta_net_faultclear(i,:) = delta_net_temp';
       %voltage_net_faultclear(i,:) = voltage_temp';
       [Results_fsolve,fval,exitflag,output]=fsolve(@(x)Fun_AEfslove_SPM(x,fault_delta_gen(i,:)',preset,system),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
       if exitflag<=0
           error('cannot find fault-clear state! \n');
       end
       delta_net_faultclear(i,:) = Results_fsolve(1:6)';
       voltage_net_faultclear(i,:) = Results_fsolve(7:12)';
       temp_ini = Results_fsolve;
    end

%% postfault
    system="postfault";
    delta0=x_fault_all(end,1:3)';
    ydc0 = x_fault_all(end,17);
    omega0=[x_fault_all(end,16)';x_fault_all(end,18)']; %1,3
    [delta_net_s,Voltage_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(temp_ini(1:6),temp_ini(7:12),delta0,preset,Basevalue,system,1e4,1e-10);
    [Results_fsolve,fval,exitflag,output]=fsolve(@(x)Fun_AEfslove_SPM(x,delta0,preset,system),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
    
    [t_postfault, x_postfault_all] = ode15s(@(t,x)f_timedomain_DAE(t,x,system),Iter.Trecover:Iter.Tunit:Iter.Ttotal,[delta0; delta_net_s; Voltage_net_s;omega0(1);ydc0;omega0(2)],options);
    %clear delta0 omega0

    out2 = [];
    for i =2: length(t_postfault)
        out = F_3M9B_SP_DAE_VFM(x_postfault_all(i,:),system);
        out2 = [out2;out(2)];
    end
    omega_VFM_postfault = out2;

%% data collection
     t_timedomain_DAE = [t_prefault(2:end);t_fault(2:end);t_postfault(2:end)];
     deltac_timedomain_DAE  = [x_prefault_all(2:end,1:3); fault_delta_gen; x_postfault_all(2:end,1:3)];
     omega_timedomain_DAE1  = [x_prefault_all(2:end,16); x_fault_all(2:end,16); x_postfault_all(2:end,16)]./omegab;
     omega_timedomain_DAE3  = [x_prefault_all(2:end,18); x_fault_all(2:end,18); x_postfault_all(2:end,18)]./omegab;
     omegac_timedomain_DAE_VFM = [omega_VFM_prefault;omega_VFM_fault;omega_VFM_postfault]./omegab;

     ydc_timedomain_DAE = [x_prefault_all(2:end,17); x_fault_all(2:end,17); x_postfault_all(2:end,17)];
     ydc_timedomain_fault = x_fault_all(2:end,17);
     ydc_timedomain_post = x_postfault_all(2:end,17);

     omegacoi_timedomain_fault = (x_fault_all(2:end,16) *preset.m(1)+ x_fault_all(2:end,18) *preset.m(3))/sum(preset.m);
     omegacoi_timedomain_post = (x_postfault_all(2:end,16) *preset.m(1)+ x_postfault_all(2:end,18) *preset.m(3))/sum(preset.m);

     omegacoi_timedomain_DAE  = (omega_timedomain_DAE1 *preset.m(1)+ omega_timedomain_DAE3 *preset.m(3))/sum(preset.m);
     omegac_timedomain_DAE1  = omega_timedomain_DAE1  - omegacoi_timedomain_DAE ;
     omegac_timedomain_DAE3  = omega_timedomain_DAE3  - omegacoi_timedomain_DAE ;

     delta_net_DAE = [x_prefault_all(2:end,4:9); fault_delta_net; x_postfault_all(2:end,4:9)];
     voltage_net_DAE = [x_prefault_all(2:end,10:15); fault_voltage_net; x_postfault_all(2:end,10:15)];

     delta_net_clear_DAE = [x_prefault_all(2:end,4:9); delta_net_faultclear; x_postfault_all(2:end,4:9)];
     voltage_net_clear_DAE = [x_prefault_all(2:end,10:15); voltage_net_faultclear; x_postfault_all(2:end,10:15)];




%% plot for DAE
% deltac
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE,deltac_timedomain_DAE(:,1),'linewidth',2,'color',[0/255 95/255 255/255]);    hold on; %blue
    plot(t_timedomain_DAE,deltac_timedomain_DAE(:,2),'linewidth',2,'color',[255/255 135/255 0/255]);    hold on; %orange
    plot(t_timedomain_DAE,deltac_timedomain_DAE(:,3),'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;  %green


    grid on;    grid minor;
    ylim([-1,2]);
    yl=ylim;
    ymin=yl(1,1);
    ymax=yl(1,2);
    xlim([20-T_before 20+T_after]);
    thetarange=[ymin,ymin,ymax,ymax];
    xl=xlim;
    xmin=xl(1,1);
    xmax=xl(1,2);
    ylab_x=xmin-(xmax-xmin)/15;
    ylab_y=(ymax+ymin)/2;
    xlab_x=(xmax+xmin)/2;
    xlab_y=ymin-(ymax-ymin)/15;
    ylim([ymin ymax]);
    trange=[Iter.Tfault,Iter.Trecover,Iter.Trecover,Iter.Tfault];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    axis([xmin xmax ymin ymax ]);
    % axis font
    ax=gca;
    ax.FontName='Arial';
    ax.FontSize=14;
    legendEntries = cell(1, 3);
    for i = 1:3
        legendEntries{i} = ['IBR' num2str(i)];
    end
    legend(legendEntries{:});
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('$\delta_c(rad)$','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);
%% omega
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);

    plot(t_timedomain_DAE,omegac_timedomain_DAE1,'linewidth',2,'color',[0/255 95/255 255/255]);    hold on;
    plot(t_timedomain_DAE,omegac_timedomain_DAE_VFM,'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
    plot(t_timedomain_DAE,omegac_timedomain_DAE3,'linewidth',2,'color',[0/255 175/255 0/255]);    hold on;
    plot(t_timedomain_DAE,omegacoi_timedomain_DAE-1,'linewidth',2,'color',[0/255 0/255 0/255]);    hold on;

    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;
    ylim([-0.02,0.03]);

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
    legendEntries = cell(1, 3);
    for i = 1:3
        legendEntries{i} = ['IBR' num2str(i)];
    end
    legend(legendEntries{:});
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('$\omega(pu)$','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);

    %% dc link
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE,ydc_timedomain_DAE,'linewidth',2,'color',[255/255 135/255 0/255]);    hold on;
    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;

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
    set(get(ax,'YLabel'),'Interpreter','latex');
    legend('DC voltage error of VFM');
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
xlabel('$\mathrm{Time\,(s)}$', 'Interpreter','latex', ...
       'FontSize',18, 'FontName','Times New Roman', ...
       'FontAngle','italic','FontWeight','bold', ...
       'Position',[xlab_x, xlab_y, -1]);

ylabel('$\Delta\!\left( v_{\mathrm{dc}} \right)^{2}\;(\mathrm{pu})$', ...
       'Interpreter','latex', 'FontSize',18, 'FontName','Times New Roman', ...
       'FontAngle','italic','FontWeight','bold', ...
       'Position',[ylab_x, ylab_y, -1]);

    
    
    %% Bus angle
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);



    plot(t_timedomain_DAE,delta_net_DAE(:,1),'linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,2),'linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,3),'linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,4),'linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,5),'linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_timedomain_DAE,delta_net_DAE(:,6),'linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;

    plot(t_fault(2:end),delta_net_faultclear(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_fault(2:end),delta_net_faultclear(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;


    
    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;
    ylim([-1,2]);
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
    legendEntries = cell(1, 6);
    for i = 1:6
        legendEntries{i} = ['Bus' num2str(i+3)];
    end
    legend(legendEntries{:});
    set(groot, 'defaultTextInterpreter','latex', ...
           'defaultAxesTickLabelInterpreter','latex', ...
           'defaultLegendInterpreter','latex');
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('$\theta_c(rad)$','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);

    %% Bus voltage
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);

    plot(t_timedomain_DAE,voltage_net_DAE(:,1),'linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,2),'linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,3),'linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,4),'linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,5),'linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_timedomain_DAE,voltage_net_DAE(:,6),'linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;

    plot(t_fault(2:end),voltage_net_faultclear(:,1),'LineStyle',':','linewidth',2,'color',[0.8500 0.3250 0.0980]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,2),'LineStyle',':','linewidth',2,'color',[0.9290 0.6940 0.1250]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,3),'LineStyle',':','linewidth',2,'color',[0.4940 0.1840 0.5560]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,4),'LineStyle',':','linewidth',2,'color',[0.4660 0.6740 0.1880]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,5),'LineStyle',':','linewidth',2,'color',[0.3010 0.7450 0.9330]);    hold on;
    plot(t_fault(2:end),voltage_net_faultclear(:,6),'LineStyle',':','linewidth',2,'color',[0.6350 0.0780 0.1840]);    hold on;

    xlim([20-T_before 20+T_after]);
    grid on;    grid minor;
    ylim([-0.1,1.2]);
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
    legendEntries = cell(1, 6);
    for i = 1:6
        legendEntries{i} = ['Bus' num2str(i+3)];
    end
    legend(legendEntries{:});
    xlabel('Time(s)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[xlab_x,xlab_y,-1]);
    ylabel('voltage(p.u.)','FontSize',18,'FontName','Times New Roman','FontAngle','italic','FontWeight','bold','position',[ylab_x,ylab_y,-1]);

%% critical energy

deltaccc = ep_set(1).xep;
delta1ccc= -preset.m(2:3)'*deltaccc(1:2)/preset.m(1);
deltaccc = [delta1ccc; deltaccc];
delta_net_cri = ep_set(1).delta_net_ep;
voltage_net_cri = ep_set(1).voltage_net_ep;
system = "postfault";
Fun_AEfslove_SPM([delta_net_cri;voltage_net_cri],deltacc,preset,system)
f_reducedstate_SPM_VFM([ep_set(1).xep;delta_net_cri;voltage_net_cri])

[Results_fsolve,fval,exitflag,output]=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc,preset,system),[delta_net_cri;voltage_net_cri],optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
%%
[Ep_tmp(1),Ep_tmp(2),Ep_tmp(3),Ep_tmp(4),Ep_tmp(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltacc',delta_net_cri,voltage_net_cri);
energy_critical = sum(Ep_tmp);

clear Ep_tmp


%% Calculate Energy
     no_duration = n_end - n_start + 1;
     Pm=preset.Pmpu;
     E=preset.Epu;
     m=preset.m;
     d=preset.d;
     Yfull_post=postfault.Yfull_mod;
     G_post=real(Yfull_post);
     B_post=imag(Yfull_post);
    % Potential energy
        [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3),Ep_tmp(4),Ep_tmp(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_timedomain_DAE(n_start,:)',delta_net_clear_DAE(n_start,:)',voltage_net_clear_DAE(n_start,:)');
        Ep_start=sum(Ep_tmp);
        Ep_loadloss_start=Ep_tmp(3);
        Ep_mag_start=Ep_tmp(2);
        Ep_mag_load_start = Ep_tmp(5);
        Ep_lossy_start = Ep_tmp(4);
        clear Ep_tmp
        Ep=zeros(no_duration,1);
        Ep_loadloss=zeros(no_duration,1);
        Ep_lossy=zeros(no_duration,1);
        Ep_mag=zeros(no_duration,1);
        Ep_mag_load=zeros(no_duration,1);
        for tm=1:size(Ep,1)
            [Ep_tmp(1),Ep_tmp(2),Ep_tmp(3),Ep_tmp(4),Ep_tmp(5)]=Fun_Cal_PotentialEnergy_SPM(preset,postfault,deltac_timedomain_DAE(n_start-1+tm,:)',delta_net_clear_DAE(n_start-1+tm,:)',voltage_net_clear_DAE(n_start-1+tm,:)');
            Ep_loadloss(tm)=Ep_tmp(3);
            Ep_mag(tm)=Ep_tmp(2);
            Ep_lossy(tm)=Ep_tmp(4);
            Ep_mag_load(tm)=Ep_tmp(5);
            Ep(tm)=sum(Ep_tmp);   clear Ep_tmp
        end

     % Kinetic energy
        Ek0=zeros(ngen,1);
        Eke=zeros(ngen,1);
        Ek0(1)=0.5*m(1)*(omegac_timedomain_DAE1(n_start)*omegab)^2;
        Ek0(2)=0.5*m(3)*(omegac_timedomain_DAE3(n_start)*omegab)^2;
        Ek0dc = 0.25*preset.C_dc*preset.kip*ydc_timedomain_DAE(n_start)^2;
        Ek_start=sum(Ek0)+Ek0dc;

    % Damping energy in iteration
        % for Potential energy error display
        Ed_non_record=zeros(no_duration,1);
        Ed_un_record=zeros(no_duration,1);
        Ep_lossy_iter=zeros(no_duration,1);
        Ep_loadloss_iter=zeros(no_duration,1);
        Ep_magload_iter=zeros(no_duration,1);
        Ep_mag_iter=zeros(no_duration,1);
        Ep_iter_record=zeros(no_duration,1);
        Ep_iter_record_withoutnet=zeros(no_duration,1);
        P_un_record=zeros(no_duration,1);
        P_non_record=zeros(no_duration,1);

        Pd_vfm_record=zeros(no_duration,1);
        Ed_vfm_iter=zeros(no_duration,1);
        Ed_coi_iter=zeros(no_duration,1);
        P_coi_record = zeros(no_duration,1);

        Ep_iter_record(1)=Ep_start;
        Ep_iter_record_withoutnet(1)=Ep_start;
        Ep_lossy_iter(1)=Ep_lossy_start;
        Ep_loadloss_iter(1) = Ep_loadloss_start;
        Ep_magload_iter(1) = Ep_mag_load_start;
        Ep_mag_iter(1)=Ep_mag_start;
        Ek=zeros(no_duration,1);
        Ek(1)=Ek_start;
        Err_step=zeros(no_duration,1);
        Err_step_exp=zeros(no_duration,1);

        Pe_pre=zeros(ngen,1);
        Pnet_pre=zeros(nbus-ngen,1);
        Qnet_pre=zeros(nbus-ngen,1);
        Pe_lossy_pre=zeros(ngen,1);  % lossy dispattive energy
        Pnet_lossy_pre=zeros(nbus-ngen,1);  % lossy dispattive energy
        Qnet_lossy_pre=zeros(nbus-ngen,1);  % lossy dispattive energy
        P_loadloss_pre=zeros(nbus-ngen,1);  % lossy dispattive energy
        Pe_mag_pre=zeros(ngen,1);   % Pe relevant to line inductance
        Qnet_mag_pre=zeros(nbus-ngen,1);   % Pe relevant to line inductance
        Pnet_mag_pre=zeros(nbus-ngen,1);   % Pe relevant to line inductance
        mag_load_P_pre=zeros(nbus-ngen,1);   % P Q load function on net
        mag_load_Q_pre=zeros(nbus-ngen,1);   % P Q load function on net
        for tm=2:no_duration
            ddelta_tmp=deltac_timedomain_DAE(n_start-1+tm,:)-deltac_timedomain_DAE(n_start-1+tm-1,:);
            dydc_tmp=ydc_timedomain_DAE(n_start-1+tm)-ydc_timedomain_DAE(n_start-1+tm-1);
            dtheta_tmp=delta_net_clear_DAE(n_start-1+tm,:)-delta_net_clear_DAE(n_start-1+tm-1,:);
            dvoltage_tmp=voltage_net_clear_DAE(n_start-1+tm,:)-voltage_net_clear_DAE(n_start-1+tm-1,:);
            Pe=zeros(ngen,1);
            Pnet=zeros(nbus-ngen,1);
            Qnet=zeros(nbus-ngen,1);
            Pe_lossy=zeros(ngen,1);  % lossy dispattive energy
            Pnet_lossy=zeros(nbus-ngen,1);  % lossy dispattive energy
            Qnet_lossy=zeros(nbus-ngen,1);  % lossy dispattive energy
            P_loadloss=zeros(nbus-ngen,1);  % lossy dispattive energy
            Pe_mag=zeros(ngen,1);   % Pe relevant to line inductance
            Qnet_mag=zeros(nbus-ngen,1);   % Pe relevant to line inductance
            Pnet_mag=zeros(nbus-ngen,1);   % Pe relevant to line inductance

            mag_load_P=zeros(nbus-ngen,1);   % P Q load function on net
            mag_load_Q=zeros(nbus-ngen,1);   % P Q load function on net

            Ed_un=zeros(ngen,1);
            Ed_un_this=zeros(ngen,1);
            Ed_non=zeros(ngen,1);
            Ed_non_this=zeros(ngen,1);
            Pd_un=zeros(ngen,1);
            Pd_non=zeros(ngen,1);


            Ep_iter_gen=zeros(ngen,1);
            Ep_iter_net=zeros((nbus-ngen),1);


            % P calculation of gen
            for i=1:ngen
                for j=1:ngen
                    ddelta=deltac_timedomain_DAE(n_start-1+tm,i)-deltac_timedomain_DAE(n_start-1+tm,j);
                    Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
                    if(i~=j)
                        Pe_lossy(i)=Pe_lossy(i)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
                        Pe_mag(i)=Pe_mag(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta);
                    end
                end
                for l=1:(nbus-ngen)
                    ddelta=deltac_timedomain_DAE(n_start-1+tm,i)-delta_net_clear_DAE(n_start-1+tm,l);
                    Pe(i)=Pe(i)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i,l+ngen)*cos(ddelta);
                    Pe_mag(i)=Pe_mag(i)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i,l+ngen)*sin(ddelta);
                    Pe_lossy(i)=Pe_lossy(i)+E(i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i,l+ngen)*cos(ddelta);
                end


            end
            %damping energy
            Ed_un(1)=preset.d(1)*omegac_timedomain_DAE1(n_start-1+tm-1,1)*omegab*ddelta_tmp(1);
            Ed_un_this(1)=preset.d(1)*omegac_timedomain_DAE1(n_start-1+tm,1)*omegab*ddelta_tmp(1);
            Pd_un(1)=-preset.d(1)*(omegac_timedomain_DAE1(n_start-1+tm-1,1)*omegab)^2;
            Ed_non(1)=preset.d(1)*(omegacoi_timedomain_DAE(n_start-1+tm-1,1)*omegab-omegab)*ddelta_tmp(1);
            Ed_non_this(1)=preset.d(1)*(omegacoi_timedomain_DAE(n_start-1+tm,1)*omegab-omegab)*ddelta_tmp(1);
            Pd_non(1)=-preset.d(1)*(omegacoi_timedomain_DAE(n_start-1+tm-1,1)*omegab-omegab)*(omegac_timedomain_DAE1(n_start-1+tm-1,1)*omegab);

            Ed_un(3)=preset.d(3)*omegac_timedomain_DAE3(n_start-1+tm-1,1)*omegab*ddelta_tmp(3);
            Ed_un_this(3)=preset.d(3)*omegac_timedomain_DAE3(n_start-1+tm,1)*omegab*ddelta_tmp(3);
            Pd_un(3)=-preset.d(3)*(omegac_timedomain_DAE3(n_start-1+tm-1,1)*omegab)^2;
            Ed_non(3)=preset.d(3)*(omegacoi_timedomain_DAE(n_start-1+tm-1,1)*omegab-omegab)*ddelta_tmp(3);
            Ed_non_this(3)=preset.d(3)*(omegacoi_timedomain_DAE(n_start-1+tm,1)*omegab-omegab)*ddelta_tmp(3);
            Pd_non(3)=-preset.d(3)*(omegacoi_timedomain_DAE(n_start-1+tm-1,1)*omegab-omegab)*(omegac_timedomain_DAE3(n_start-1+tm-1,1)*omegab);



            
            clear i j l
            % P calculation of Bus
            for i=1:(nbus-ngen)
                for j=1:ngen
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-deltac_timedomain_DAE(n_start-1+tm,j);
                    Pnet(i)=Pnet(i)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*B_post(i+ngen,j)*sin(ddelta)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*G_post(i+ngen,j)*cos(ddelta);
                    Pnet_lossy(i)=Pnet_lossy(i)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*G_post(i+ngen,j)*cos(ddelta);
                    Pnet_mag(i)=Pnet_mag(i)+voltage_net_clear_DAE(n_start-1+tm,i)*E(j)*B_post(i+ngen,j)*sin(ddelta);
                end
                for l=1:(nbus-ngen)
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-delta_net_clear_DAE(n_start-1+tm,l);
                    Pnet(i)=Pnet(i)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*sin(ddelta)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*cos(ddelta);
                    if(i~=l)
                        Pnet_lossy(i)=Pnet_lossy(i)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*cos(ddelta);
                        Pnet_mag(i)=Pnet_mag(i)+voltage_net_clear_DAE(n_start-1+tm,i)*voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*sin(ddelta);
                    else
                        P_loadloss(i) = voltage_net_clear_DAE(n_start-1+tm,i)^2*G_post(i+ngen,i+ngen);
                    end
                end
                for h=1:size(preset.Sload,1)
                    if (preset.Sload(h,1)==postfault.Transform(i+ngen))
                         Pnet(i)=Pnet(i)+preset.Sload(h,2);
                         mag_load_P(i)=mag_load_P(i)+preset.Sload(h,2);
                    end
                end
            end
            clear i j l h
            % Q/V calculation of Bus
            for i=1:(nbus-ngen)
                for j=1:ngen
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-deltac_timedomain_DAE(n_start-1+tm,j);
                    Qnet(i)=Qnet(i)-E(j)*B_post(i+ngen,j)*cos(ddelta)+E(j)*G_post(i+ngen,j)*sin(ddelta);
                    Qnet_mag(i)=Qnet_mag(i)-E(j)*B_post(i+ngen,j)*cos(ddelta);
                    Qnet_lossy(i)=Qnet_lossy(i)+E(j)*G_post(i+ngen,j)*sin(ddelta);
                end
                for l=1:(nbus-ngen)
                    ddelta=delta_net_clear_DAE(n_start-1+tm,i)-delta_net_clear_DAE(n_start-1+tm,l);
                    Qnet(i)=Qnet(i)-voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*cos(ddelta)+voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*sin(ddelta);
                    Qnet_mag(i)=Qnet_mag(i)-voltage_net_clear_DAE(n_start-1+tm,l)*B_post(i+ngen,l+ngen)*cos(ddelta);
                    Qnet_lossy(i)=Qnet_lossy(i)+voltage_net_clear_DAE(n_start-1+tm,l)*G_post(i+ngen,l+ngen)*sin(ddelta);
                end
                for h=1:size(preset.Sload,1)
                    if (preset.Sload(h,1)==postfault.Transform(i+ngen))
                           Qnet(i)=Qnet(i)+preset.Sload(h,3)/voltage_net_clear_DAE(n_start-1+tm,i);
                           mag_load_Q(i)=mag_load_Q(i)+preset.Sload(h,3)/voltage_net_clear_DAE(n_start-1+tm,i);
                    end
                end
            end  
            clear i j l h

            if tm==2
                Pe_pre=Pe;
                Pnet_pre=Pnet;
                Qnet_pre=Qnet;
                Pe_lossy_pre=Pe_lossy;  % lossy dispattive energy
                Pnet_lossy_pre=Pnet_lossy;  % lossy dispattive energy
                Qnet_lossy_pre=Qnet_lossy;  % lossy dispattive energy
                P_loadloss_pre=P_loadloss;  % lossy dispattive energy
                Pe_mag_pre=Pe_mag;   % Pe relevant to line inductance
                Qnet_mag_pre=Qnet_mag;   % Pe relevant to line inductance
                Pnet_mag_pre=Pnet_mag;   % Pe relevant to line inductance
                mag_load_P_pre=mag_load_P;
                mag_load_Q_pre=mag_load_Q;
            end


            for i=1:ngen
                Ep_iter_gen(i)=(Pm(i)-Pe(i)+Pm(i)-Pe_pre(i))/2*ddelta_tmp(i);
            end
            Ek(tm)=Ek(tm)+0.5*m(1)*(omegac_timedomain_DAE1(n_start+tm-1)*omegab)^2;
            Ek(tm)=Ek(tm)+0.5*m(3)*(omegac_timedomain_DAE3(n_start+tm-1)*omegab)^2;
            Ekdc = 0.25*preset.C_dc*preset.kip*ydc_timedomain_DAE(n_start+tm-1)^2;
            Ek(tm)= Ek(tm) + Ekdc;
           
            for i=1:nbus-ngen
                Ep_iter_net(i)=(Pnet(i)+Pnet_pre(i))/2*dtheta_tmp(i) + (Qnet(i)+Qnet_pre(i))/2*dvoltage_tmp(i);
            end
            if tm==1000
                Pm(2)-Pe(2)
            end

            %damping energy on VFM
            Pd_vfm_record(tm) = -preset.kpp*(Pm(2)-Pe(2))^2;
            Ed_vfm_iter(tm) = Ed_vfm_iter(tm-1) + preset.kpp*((Pm(2)-Pe(2))+ (Pm(2)-Pe_pre(2)))/2*preset.C_dc/2*dydc_tmp;

            %damping brought by coi
            Ed_coi_iter(tm) = Ed_coi_iter(tm-1) - preset.C_dc/2* ( (omegacoi_timedomain_DAE(n_start-1+tm,1)+omegacoi_timedomain_DAE(n_start-1+tm-1,1)-2)/2*omegab ) * dydc_tmp;
            P_coi_record(tm) = (Pm(2)-Pe(2))*(omegacoi_timedomain_DAE(n_start-1+tm,1)*omegab-omegab);




            Ep_lossy_iter(tm)=Ep_lossy_iter(tm-1)+ddelta_tmp*(Pe_lossy+Pe_lossy_pre)/2+dtheta_tmp*(Pnet_lossy+Pnet_lossy_pre)/2+dvoltage_tmp*(Qnet_lossy+Qnet_lossy_pre)/2;   %losses energy on network exclude self-susceptance
            Ep_mag_iter(tm)=Ep_mag_iter(tm-1)+ddelta_tmp*(Pe_mag+Pe_mag_pre)/2+dtheta_tmp*(Pnet_mag+Pnet_mag_pre)/2+dvoltage_tmp*(Qnet_mag+Qnet_mag_pre)/2;% %magnatic potential energy
            Ep_loadloss_iter(tm)=Ep_loadloss_iter(tm-1)+dtheta_tmp*(P_loadloss+P_loadloss_pre)/2; %R-type load losses containing self-susceptance
            
            Ep_iter_record(tm)=Ep_iter_record(tm-1)-sum(Ep_iter_gen) + sum(Ep_iter_net); %sum Ep Iteration
            Ep_iter_record_withoutnet(tm)= Ep_iter_record(tm-1)-sum(Ep_iter_gen); %sum Ep Iteration

            Ep_magload_iter(tm) = Ep_magload_iter(tm-1)+dtheta_tmp*(mag_load_P+mag_load_P_pre)/2+dvoltage_tmp*(mag_load_Q_pre+mag_load_Q)/2;



            Ed_non_record(tm)=Ed_non_record(tm-1)+(sum(Ed_non)+sum(Ed_non_this))/2;
            Ed_un_record(tm)=Ed_un_record(tm-1)+(sum(Ed_un)+sum(Ed_un_this))/2;

            Err_step(tm)=Ep_iter_record(tm)+Ed_un_record(tm)+Ed_non_record(tm)+Ek(tm)-Ek(1)+Ed_vfm_iter(tm)+Ed_coi_iter(tm); %sum energy
            Err_step_exp(tm)=Ep(tm)+Ed_un_record(tm)+Ed_non_record(tm)+Ek(tm)-Ek(1)+Ed_vfm_iter(tm)+Ed_coi_iter(tm); %sum energy

            P_un_record(tm)=sum(Pd_un);
            P_non_record(tm)=sum(Pd_non);








            % Trapezoidal Int
            Pe_pre=Pe;
            Pnet_pre=Pnet;
            Qnet_pre=Qnet;
            Pe_lossy_pre=Pe_lossy;  % lossy dispattive energy
            Pnet_lossy_pre=Pnet_lossy;  % lossy dispattive energy
            Qnet_lossy_pre=Qnet_lossy;  % lossy dispattive energy
            P_loadloss_pre=P_loadloss;  % lossy dispattive energy
            Pe_mag_pre=Pe_mag;   % Pe relevant to line inductance
            Qnet_mag_pre=Qnet_mag;   % Pe relevant to line inductance
            Pnet_mag_pre=Pnet_mag;   % Pe relevant to line inductance
            mag_load_P_pre=mag_load_P;
            mag_load_Q_pre=mag_load_Q;
        end                                                  

    %% plot

    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_loadloss,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_loadloss_iter,'color','g','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential Ray','Iter');
    title('Resistive lossy energy on load');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_lossy,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_lossy_iter,'color','g','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Calculation','Iter');
    title('lossy energy in grid');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep_mag,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_mag_iter,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Potential','Iter');
    title('magnetic energy');
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end), Ep_loadloss-Ep_loadloss_iter,'color','g','LineWidth',2);  hold on; 
    plot(t_timedomain_DAE(n_start:n_end), Ep_lossy-Ep_lossy_iter,'color','y','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end), Ep-Ep_iter_record,'color','b','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end), Ep_loadloss-Ep_loadloss_iter+Ep_lossy-Ep_lossy_iter,'color','r','LineWidth',2,'LineStyle','--');  hold on;
    plot(t_timedomain_DAE(n_start:n_end), Ep_iter_record-Ep_iter_record_withoutnet,'color','k','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('load Lossy part err','network Lossy part err','Total err','Total lossy err','network potential error');
    title('Potential Energy Error');  

    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ed_un_record,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ed_vfm_iter,'color','b','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ed_coi_iter,'color','g','LineWidth',2,'LineStyle','-');  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ed_coi_iter+Ed_vfm_iter+Ed_un_record,'color','k','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('VSG damping','vfm damping','COI damping','total');
    title('Damping Energy');


%%
    figure;
    grid on;
    hold on;
    
    set(gca,'FontSize',16); 
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Ep+Ek,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep_iter_record+Ek,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Ep+Ek+Ed_coi_iter,'color','k','LineWidth',2,'LineStyle','-');  hold on;
    plot([t_timedomain_DAE(n_start);t_timedomain_DAE(end)],[2.93;2.93],'k','LineWidth',2,'LineStyle','-'); hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    axis([19 24 -0.1 4.5]);
    legend('Expression','Iteration','Expression + COI damping');
    title('Energy Function'); 
%%


    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),P_un_record,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Pd_vfm_record,'color','b','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),P_coi_record,'color','g','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),P_un_record+Pd_vfm_record+P_coi_record,'color','k','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('VSG damping','vfm damping','COI damping','total');
    title('Derivative of damping energy'); 
%%

    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    plot(t_timedomain_DAE(n_start:n_end),Err_step_exp,'color','r','LineWidth',2);  hold on;
    plot(t_timedomain_DAE(n_start:n_end),Err_step,'color','b','LineWidth',2,'LineStyle','--');  hold on;
    ax = gca; yl=ax.YLim; ymin=yl(1,1); ymax=yl(1,2); thetarange=[ymin,ymin,ymax,ymax];
    fill(trange,thetarange,[.9805 .7031 .6797], 'linestyle', 'none', 'FaceAlpha',0.5);
    legend('Expression','Iteration');
    title('Total Energy Change');  

%%
    figure;
    set(gca,'position',[0.115,0.12,0.815,0.84]);
    set(gcf,'position',[60 200 600 450]);
    grid on; hold on;
    plot(ydc_timedomain_fault,omegacoi_timedomain_fault./omegab-1,'color','r','LineWidth',2,'LineStyle','--');  hold on;
    plot(ydc_timedomain_post,omegacoi_timedomain_post./omegab-1,'color','b','LineWidth',2,'LineStyle','--');  hold on;


%% ode function
function dfdt = f_timedomain(t,x)
    dfdt = F_3M9B_MR_ODE(x);
end

function dfdt = f_timedomain_DAE(t,x,system)
    dfdt = F_3M9B_SP_DAE_VFM(x,system);
end
function dfdt = f_timedomain_DAE_withflove(t,x,net_guess,system)
    net_value=fsolve(@(x)Fun_AEfslove_SPM(x,delta0,preset,system),net_guess,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',1e-9));
    dfdt = F_3M9B_SP_ODE(x,net_value,system);
end





