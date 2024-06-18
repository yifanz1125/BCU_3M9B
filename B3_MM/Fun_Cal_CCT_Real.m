%% function: calculate real CCT along fault-on trajectory
%% unstable criterion: 
function [CCT,Exit_thetac,Exit_omegac,Exit_theta,Exit_omega,flag_CCT,Traj_Stb,Traj_Unstb]=Fun_Cal_CCT_Real(fault,postfault,preset,Basevalue,CCT_ref)
    cycle=size(fault.traj.omega,1);
    ngen=size(postfault.Yred,1);
    omegab=Basevalue.omegab;
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    flag_CCT=0;
    Tunitmin=fault.traj.Tunit;
    Tfaultlen=fault.traj.Tlength;
    thetac_fault=fault.traj.thetac;
    omega_fault=fault.traj.omega;
    
%% Change unit settings
    Changeunit.n_stepchange=5; % numbers of units
    Changeunit.candidate=[1e3;1e2;10;5;1];
    if(Tfaultlen<=2*Tunitmin*Changeunit.candidate(1))
        error('Unsuitable Tunit selected, choose smaller T_unit or regenerate fault traj with longer time');
    end
%% Stable & Unstable settings
    Stb.n_current=fix(CCT_ref/Tunitmin);    % start from reference CCT
    Stb.n_stbmax=1;
    Stb.n_unstbmin=cycle;
%% Postfault traj iteration settings
    Tunit=1e-4;
    Tpostmax=200;    % maximum time for postfault iteration
%% Prejudge the CCT_ref
    delta0=thetac_fault(Stb.n_current,:);
    omega0=omega_fault(Stb.n_current,:);
    [theta_RK4,omega_RK4,thetac_RK4,omegacoi_RK4,Pe,cycle,flag_unstb]=Fun_TrajIter_StableCheck_SRF(Tpostmax,Tunit,postfault,preset,delta0,omega0,omegab);
    if(flag_unstb==0)
        Stb.n_stbmax=Stb.n_current;
        Stb.n_current=(fix(Stb.n_stbmax/Changeunit.candidate(1))+1)*Changeunit.candidate(1);
        Traj_Stb.thetac=fault.traj.thetac;
        Traj_Stb.omegac=fault.traj.omegac;
        Traj_Stb.theta=fault.traj.theta;
        Traj_Stb.omega=fault.traj.omega;
        
    else
        Stb.n_unstbmin=Stb.n_current;
        Stb.n_current=(fix(Stb.n_unstbmin/Changeunit.candidate(1)))*Changeunit.candidate(1);
    end
%% Change unit
    for n_step=1:Changeunit.n_stepchange
        Changeunit.Currentunit=Changeunit.candidate(n_step);
        flag_changeunit=0;
    %% Current unit scan
        while(flag_changeunit==0)
            delta0=thetac_fault(Stb.n_current,:);
            omega0=omega_fault(Stb.n_current,:);
            [theta_RK4,omega_RK4,thetac_RK4,omegacoi_RK4,Pe,cycle,flag_unstb]=Fun_TrajIter_StableCheck_SRF(Tpostmax,Tunit,postfault,preset,delta0,omega0,omegab);
            
            if(flag_unstb==0)
                if(Stb.n_stbmax<Stb.n_current)
                    Stb.n_stbmax=Stb.n_current;
                    Traj_Stb.thetac=thetac_RK4(1:cycle,:);
                    Traj_Stb.omegac=omega_RK4-omegacoi_RK4*ones(1,ngen);
                    Traj_Stb.theta=theta_RK4(1:cycle,:);
                    Traj_Stb.omega=omega_RK4(1:cycle,:);
                end
            else
                if(Stb.n_unstbmin>Stb.n_current)    % for initialization
                    Stb.n_unstbmin=Stb.n_current;                    
                    Traj_Unstb.thetac=thetac_RK4(1:cycle,:);
                    Traj_Unstb.omegac=omega_RK4(1:cycle,:)-omegacoi_RK4(1:cycle)*ones(1,ngen);
                    Traj_Unstb.theta=theta_RK4(1:cycle,:);
                    Traj_Unstb.omega=omega_RK4(1:cycle,:);
                end
            end
            if(Stb.n_unstbmin-Stb.n_stbmax<=Changeunit.Currentunit)
                flag_changeunit=1;
                if(n_step<Changeunit.n_stepchange)
                    Changeunit.Nextunit=Changeunit.candidate(n_step+1);
                    Stb.n_current=(fix(Stb.n_stbmax/Changeunit.Nextunit)+1)*Changeunit.Nextunit;
                else
                    CCT=Stb.n_stbmax*Changeunit.candidate(n_step)*Tunitmin;
                    flag_CCT=1;
                end
            end
            if(Stb.n_unstbmin<=Stb.n_stbmax)
                error('Stb.n_unstbmin<=Stb.n_stbmax');
            end
            if(flag_changeunit==0)
                if(flag_unstb==0)
                    Stb.n_current=Stb.n_current+Changeunit.Currentunit;
                else
                    Stb.n_current=Stb.n_current-Changeunit.Currentunit;
                if(Stb.n_current>cycle)
                    error('Stb.n_current>cycle');
                elseif(Stb.n_current<1)
                    error('Stb.n_current<1');
                end
            end
        end
    end
    
    end
        if(flag_CCT==1)
        Exit_thetac=fault.traj.thetac(Stb.n_stbmax,:);
        Exit_omegac=fault.traj.omegac(Stb.n_stbmax,:);
        Exit_theta=fault.traj.theta(Stb.n_stbmax,:);
        Exit_omega=fault.traj.omega(Stb.n_stbmax,:);
    else
        error('No CCT found!');
    end
end