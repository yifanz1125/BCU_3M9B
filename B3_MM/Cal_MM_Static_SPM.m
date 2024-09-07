%% Input: preset (m,d,Pm,Epu,xe1,flagxd(0-default)) fault (faultline faultposition) Case (data for matpower)
%% Output: pfdata--powerflow data; netdata--Yorg(admittance matrix without load); prefault/fault/postfault (Yred,SEP)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%% Energy estimation
    preset.PathEnergyCal=10; % 0: Ray approximation  n: N-segment trap approximation (-1)--neglect this part
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
    preset.PloadZIP = [0 0 1]; % Z I P
    preset.QloadZIP = [0 0 1]; % Z I P
% 9 bus sys===========================================
%%%%%%% 思考一下怎么把这几项融到powerflow计算中 %%%%%%%%%
    preset.Pmpu=[0.8980;1.3432;0.9419];
    preset.xd1=[0.0608;0.1198;0.1813];
    preset.Epu=[1.1083;1.1071;1.0606];

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
    preset.faultline=[9;6];  % [Frombus Tobus]
    preset.faultposition=0;
    fault.faultline=preset.faultline;
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
    Yload=zeros(pfdata.bus.numload,5);
    Iload=zeros(pfdata.bus.numload,3);
    Sload=zeros(pfdata.bus.numload,3);
    Yload(:,1)=pfdata.load.no;
    Iload(:,1)=pfdata.load.no;
    Sload(:,1)=pfdata.load.no;
    Y=netdata.Y_org;
    Y_forR=netdata.Y_org;
    for i=1:pfdata.bus.numload
        PLpu=pfdata.load.PQ(i,1)/Basevalue.Sbase;
        QLpu=pfdata.load.PQ(i,2)/Basevalue.Sbase;
        VLpu=pfdata.load.voltage(i,1);
        Yload(i,2)= PLpu*preset.PloadZIP(1)/(VLpu^2);
        Yload(i,3)=-1*QLpu*preset.QloadZIP(1)/(VLpu^2);
        Yload(i,4)= PLpu/(VLpu^2);
        Yload(i,5)=-1*QLpu/(VLpu^2);
        Iload(i,2)= PLpu*preset.PloadZIP(2)/(VLpu);
        Iload(i,3)=-1*QLpu*preset.QloadZIP(2)/(VLpu);
        Sload(i,2)= PLpu*preset.PloadZIP(3);
        Sload(i,3)= QLpu*preset.QloadZIP(3);
    end
    clear PLpu QLpu VLpu
    for i=1:pfdata.bus.numload	% search for load bus
        Loadbus=Yload(i,1);
        for j=1:pfdata.bus.num  %   scan all bus
            if(j==Loadbus)
                Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
                Y_forR(j,j)=Y_forR(j,j)+complex(Yload(i,4),Yload(i,5));
            end
        end
    end
    prefault.Yfull=Y;
    prefault.Yfull_forR=Y_forR;
    preset.Yload = Yload;
    preset.Iload = Iload;
    preset.Sload = Sload;

    preset.genno = pfdata.gen.no;
    
    clear Y Loadbus Y_forR
    % Reduced Network Admittance of Prefault
    [prefault.Yred,prefault.Ynn,prefault.Ynr,prefault.Yrn,prefault.Yrr]=Fun_Yfull2Yred(prefault.Yfull_forR,pfdata,0);
    [prefault.Yfull_mod,prefault.Transform] = Fun_Yfull2Yfull(prefault.Yfull,pfdata,0);
    %% Structure Preserved Admittance of Fault
        Yfull_pre=prefault.Yfull;
        fault.Yfull=Yfull_pre;
        fault.Yfull_mod2 = fault.Yfull;
        fault.Yfull_mod2(fault.faultbus,fault.faultbus)=fault.Yfull_mod2(fault.faultbus,fault.faultbus)+1e6;
        fault.Yfull(:,fault.faultbus)=[];
        fault.Yfull(fault.faultbus,:)=[];

        fault.Yfull_forR=prefault.Yfull_forR;
        fault.Yfull_forR(:,fault.faultbus)=[];
        fault.Yfull_forR(fault.faultbus,:)=[];
        [fault.Yred,fault.Ynn,fault.Ynr,fault.Yrn,fault.Yrr]=Fun_Yfull2Yred(fault.Yfull_forR,pfdata,[1,fault.faultbus]);
        [fault.Yfull_mod,fault.Transform] = Fun_Yfull2Yfull(fault.Yfull_forR,pfdata,[1,fault.faultbus]); % during fault pure impedance %Fun_Yfull2Yfull(fault.Yfull,pfdata,[1,fault.faultbus]);
        [fault.Yfull_mod2,fault.Transform2] = Fun_Yfull2Yfull(fault.Yfull_mod2,pfdata,0);
    %% Structure Preserved Admittance of Postfault
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
        Y_forR=postfault.Yorg;
        Y=postfault.Yorg;
        for i=1:pfdata.bus.numload	% search for load bus
            Loadbus=Yload(i,1);
            for j=1:pfdata.bus.num  %   scan all bus
                if(j==Loadbus)
                    Y(j,j)=Y(j,j)+complex(Yload(i,2),Yload(i,3));
                    Y_forR(j,j)=Y_forR(j,j)+complex(Yload(i,4),Yload(i,5));
                end
            end
        end
        postfault.Yfull=Y;
        postfault.Yfull_forR=Y_forR;
        clear Y Loadbus Y_forR
        [postfault.Yfull_mod,postfault.Transform] = Fun_Yfull2Yfull(postfault.Yfull,pfdata,0);
        [postfault.Yred,postfault.Ynn,postfault.Ynr,postfault.Yrn,postfault.Yrr]=Fun_Yfull2Yred(postfault.Yfull_forR,pfdata,0);
        nbus=size(postfault.Yfull,1);
%% SEP calculation
    %% Adopt newton method
    if(preset.EquCal==1)
%         delta0=zeros(ngen,1);
%         [prefault.SEP_delta,prefault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(prefault.Yfull,preset.Pmpu,preset.Epu,preset.m,preset.d,delta0,0,Basevalue.omegab,1e4,1e-8);
%         if(flag_iter~=1)
%             error(' prefault SEP iteration cannot converage!\n');
%         end
%         clear flag_iter n_iter delta0
%         [postfault.SEP_delta,postfault.SEP_omegapu,flag_iter,n_iter]=Fun_SEPiteration(postfault.Yfull,preset.Pmpu,preset.Epu,preset.m,preset.d,prefault.SEP_delta,prefault.SEP_omegapu-1,Basevalue.omegab,1e4,1e-8);
%         if(flag_iter~=1)
%             error(' postfault SEP iteration cannot converage!\n');
%         end
%         clear flag_iter n_iter
    %% Adopt fsolve
    elseif(preset.EquCal==2)
        x_init=[zeros(ngen,1); zeros((nbus-ngen),1); ones((nbus-ngen),1)];
        Results_fsolve=fsolve(@(x)Fun_SEPfslove_SPM(x,preset,prefault,Basevalue),x_init,options);
        prefault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        deltacoi=delta_tmp'*preset.m/sum(preset.m);
        prefault.SEP_delta=delta_tmp-deltacoi;
        prefault.net_delta=Results_fsolve((ngen+1):nbus)-deltacoi;
        prefault.net_voltage=Results_fsolve((nbus+1):(2*nbus-ngen));
        clear delta_tmp Results_fsolve
        x_init(ngen)=(prefault.SEP_omegapu-1)*Basevalue.omegab;
        for i=1:ngen-1
            x_init(i)=prefault.SEP_delta(i)-prefault.SEP_delta(ngen);
        end
        for i=1:(nbus-ngen)
            x_init(i+ngen)=prefault.net_delta(i)-prefault.SEP_delta(ngen);
        end
        for i=1:(nbus-ngen)
            x_init(i+nbus)=prefault.net_voltage(i);
        end
        Results_fsolve=fsolve(@(x)Fun_SEPfslove_SPM(x,preset,postfault,Basevalue),x_init,options);
        postfault.SEP_omegapu=Results_fsolve(ngen)/Basevalue.omegab+1;
        delta_tmp=[Results_fsolve(1:ngen-1);0];
        deltacoi=delta_tmp'*preset.m/sum(preset.m);
        postfault.SEP_delta=delta_tmp-deltacoi;
        postfault.net_delta=Results_fsolve((ngen+1):nbus)-deltacoi;
        postfault.net_voltage=Results_fsolve((nbus+1):(2*nbus-ngen));
        clear delta_tmp Results_fsolve
    end
%     [prefault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(prefault,preset,prefault.SEP_delta,(prefault.SEP_omegapu-1)*Basevalue.omegab);
%     if(flag_SEPerr==1)
%         error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
%     end
%     [postfault.SEP_Perr,flag_SEPerr]=Fun_SEPcheck(postfault,preset,postfault.SEP_delta,(postfault.SEP_omegapu-1)*Basevalue.omegab);
%     if(flag_SEPerr==1)
%         error('SEP calculation error: the SEP solution cannot equilibrium condition! \n')
%     end
    clear flag_SEPerr
    
    clear i j