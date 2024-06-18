%% Input: preset (m,d,Pm,Epu,xe1,flagxd(0-default)) fault (faultline faultposition) Case (data for matpower)
%% Output: pfdata--powerflow data; netdata--Yorg(admittance matrix without load); prefault/fault/postfault (Yred,SEP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%% Energy estimation
    preset.PathEnergyCal=0; % 0: Ray approximation  n: N-segment trap approximation (-1)--neglect this part
%% Equilibrium calculation method
    preset.EquCal=2;    % 1--Newton method 2--fsolve method
%% fsolve preset
    options = optimset('TolFun',1e-50,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',1e-9);
    options.StepTolerance = 1e-50;
%% Generator parameters
% 9 bus sys===========================================
    Basevalue.omegab=2*pi*60;
    preset.m=[0.1254;0.034;0.016];  %[10;10;10]/omegab;
    preset.d=preset.m.*[0.2;0.2;0.2];
% % 39 bus sys=========================================
%     Basevalue.omegab=2*pi*50;
%     preset.m=[42.0;30.3;35.8;28.6;26;34.8;26.4;24.3;34.5;50]*2/(Basevalue.omegab);
%     preset.d=preset.m*0.05;
%     no_GFM=[3;4;10];
%     for i=1:size(no_GFM,1)
%         preset.d(no_GFM(i))=preset.m(no_GFM(i))*0.8;
% %     preset.d(no_GFM(i))=50/Basevalue.omegab;
%     end
% %     preset.d(2)=preset.m(2)*0.7;
% %     preset.d(2)=50/Basevalue.omegab;
%     clear no_GFM
%     
%     DT=sum(preset.d,1);
%     HT=sum(preset.m,1);
%%%%%%% 思考一下怎么把这几项融到powerflow计算中 %%%%%%%%%
% 9 bus sys===========================================
    preset.Pmpu=[0.8980;1.3432;0.9419];
    preset.xd1=[0.0608;0.1198;0.1813];
    preset.Epu=[1.1083;1.1071;1.0606];

% 39 bus sys=========================================
% order: bus no.30 to bus no.39
%     preset.Pmpu=[2.50;6.77871;6.5;6.32;5.08;6.5;5.6;5.4;8.3;10];
%     preset.xd1=[0.031;0.0697;0.0531;0.0436;0.132;0.05;0.049;0.057;0.057;0.006];
%     preset.Epu=[1.0929;1.1966;1.1491;1.0808;1.3971;1.1910;1.1394;1.0709;1.1368;1.0368];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ngen=size(preset.m,1);
    DHri=roundn(preset.d./preset.m,-3);
    flag_uniform=1;
    for i=2:ngen
        if(preset.d(i)/preset.m(i)~=DHri(1))
            flag_uniform=0;
        end
    end
    preset.flag_uniform=flag_uniform;
    clear flag_uniform DHri
%% Powerflow parameters
    path_matdata='C:\Users\yz7521\OneDrive - Imperial College London\BCU Code\BCU_3M9B\C1_Matpower\matpower7.1\data';
    addpath(genpath(path_matdata));
    Case=case9_v2;
%     Case=case39_modified;
    Basevalue.Sbase=Case.baseMVA;
%% run matpower--Matpower powerflow ignore the damping effect, and the results are used to calculate equivalent load
    path_matpower='C:\Users\yz7521\OneDrive - Imperial College London\BCU Code\BCU_3M9B\C1_Matpower\matpower7.1';
    addpath(genpath(path_matpower));
    pfdata=Fun_ResultBack(Case);
    if(ngen~=pfdata.bus.numgen)
        error('Powerflow data and generators data not match');
    end
    preset.flagxd=0;    % 0--consider xd' in network already
    % internal EMF calculation
    EMF=Fun_Cal_GenEMF(preset.flagxd,pfdata,preset.xd1);    
    clear path_matpower path_matdata
%% Fault settings
%     preset.flagfault=1; % 1--fault occur
    preset.faultline=[9;6];  % [Frombus Tobus]
    preset.faultposition=0;
    fault.faultline=preset.faultline;
%     fault.flagfault=preset.flagfault; % 1--fault occur
    fault.faultbus=preset.faultline(preset.faultposition+1,1);
%% Network 
    %% Add Reactance Xd' into Network
    if(preset.flagxd==0)
        pfdata.branch.RXB_xd=pfdata.branch.RXB;
        fprintf('WARNING: xd already included in casefile!\n');
    else    % add xd' into admittance matrix
        pfdata.branch.RXB_xd=pfdata.branch.RXB;
        for i=1:pfdata.bus.numgen
        no_gen=pfdata.gen.no(i,1);
            for j=1:size(pfdata.branch.RXB_xd,1)
                if(pfdata.branch.RXB_xd(j,1)==no_gen||pfdata.branch.RXB_xd(j,2)==no_gen)
                    pfdata.branch.RXB_xd(j,4)=pfdata.branch.RXB_xd(j,4)+preset.xd1(i);
                end
            end
        end
    end
    %% Transfer RXB into Structure-preserved admittance matrix
        netdata.Y_org=Fun_RXB2Yfull(pfdata.branch.RXB_xd,pfdata);   % admittance without load
    %% Add Passive load into Network
    Yload=zeros(pfdata.bus.numload,3);
    Yload(:,1)=pfdata.load.no;
    Y=netdata.Y_org;
    for i=1:pfdata.bus.numload
        PLpu=pfdata.load.PQ(i,1)/Basevalue.Sbase;
        QLpu=pfdata.load.PQ(i,2)/Basevalue.Sbase;
        VLpu=pfdata.load.voltage(i,1);
        Yload(i,2)=PLpu/(VLpu^2);
        Yload(i,3)=-1*QLpu/(VLpu^2);
    end
    clear PLpu QLpu VLpu
    for i=1:pfdata.bus.numload	% search for load bus
        Loadbus=Yload(i,1);
        for j=1:pfdata.bus.num  %   scan all bus
            if(j==Loadbus)
                Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
            end
        end
    end
    prefault.Yfull=Y;
    
    clear Y Loadbus
    %% Reduced Network Admittance of Prefault
        prefault.Yred=Fun_Yfull2Yred(prefault.Yfull,pfdata,0);
    %% Reduced Network Admittance of Fault
        Yfull_pre=prefault.Yfull;
        fault.Yfull=Yfull_pre;
        fault.Yfull(:,fault.faultbus)=[];
        fault.Yfull(fault.faultbus,:)=[];
        fault.Yred=Fun_Yfull2Yred(fault.Yfull,pfdata,[1,fault.faultbus]);
    %% Reduced Network Admittance of Postfault
        pfdata.branch.RXB_postfault=pfdata.branch.RXB_xd;
        for i=1:size(pfdata.branch.RXB_postfault,1)
            no1=pfdata.branch.RXB_postfault(i,1);
            no2=pfdata.branch.RXB_postfault(i,2);
            if((fault.faultline(1)==no1&&fault.faultline(2)==no2)||(fault.faultline(1)==no2&&fault.faultline(2)==no1))
                pfdata.branch.RXB_postfault(i,:)=[];
                break;
            end
        end
        clear no1 no2
        postfault.Yorg=Fun_RXB2Yfull(pfdata.branch.RXB_postfault,pfdata);
        Y=postfault.Yorg;
        for i=1:pfdata.bus.numload	% search for load bus
        Loadbus=Yload(i,1);
        for j=1:pfdata.bus.num  %   scan all bus
            if(j==Loadbus)
                Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
            end
        end
        end
        postfault.Yfull=Y;
        clear Y Loadbus
        postfault.Yred=Fun_Yfull2Yred(postfault.Yfull,pfdata,0);
%% SEP calculation
    %% Adopt newton method
    if(preset.EquCal==1)
        delta0=zeros(ngen,1);
        [prefault.SEP_delta,prefault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(prefault.Yred,preset.Pmpu,preset.Epu,preset.m,preset.d,delta0,0,Basevalue.omegab,1e4,1e-8);
        if(flag_iter~=1)
            error(' prefault SEP iteration cannot converage!\n');
        end
        clear flag_iter n_iter delta0
        [postfault.SEP_delta,postfault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(postfault.Yred,preset.Pmpu,preset.Epu,preset.m,preset.d,prefault.SEP_delta,prefault.SEP_omegapu-1,Basevalue.omegab,1e4,1e-8);
        if(flag_iter~=1)
            error(' postfault SEP iteration cannot converage!\n');
        end
        clear flag_iter n_iter
    %% Adopt fsolve
    elseif(preset.EquCal==2)
        deltaomega_init=zeros(ngen,1);
        Results_fsolve=fsolve(@(delta_omega)Fun_SEPfslove(delta_omega,preset,prefault,Basevalue),deltaomega_init,options);
        prefault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        prefault.SEP_delta=delta_tmp-delta_tmp'*preset.m/sum(preset.m);
        clear delta_tmp Results_fsolve
        deltaomega_init(ngen)=(prefault.SEP_omegapu-1)*Basevalue.omegab;
        for i=1:ngen-1
            deltaomega_init(i)=prefault.SEP_delta(i)-prefault.SEP_delta(ngen);
        end
        Results_fsolve=fsolve(@(delta_omega)Fun_SEPfslove(delta_omega,preset,postfault,Basevalue),deltaomega_init,options);
        postfault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        postfault.SEP_delta=delta_tmp-delta_tmp'*preset.m/sum(preset.m);
        clear delta_tmp Results_fsolve deltaomega_init
    end
    [prefault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(prefault,preset,prefault.SEP_delta,(prefault.SEP_omegapu-1)*Basevalue.omegab);
    if(flag_SEPerr==1)
        error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
    end
    [postfault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,postfault.SEP_delta,(postfault.SEP_omegapu-1)*Basevalue.omegab);
    if(flag_SEPerr==1)
        error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
    end
    clear flag_SEPerr
    
    clear i j