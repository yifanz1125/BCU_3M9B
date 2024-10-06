deltac_escape=escape.deltac;
theta_escape=escape.theta;
voltage_escape=escape.voltage;
%% Settings
    Tunit=1e-4; % time unit for iteration
    n_itermax=20;    % maximum steps in one iteration procedure
    norm_Tol=1e-5;  % Tolerance set for MGP identification
    n_MGPtraj=0;    % counter for trajectory numbers
    n_MGPtrajmax=1000;  % maximum trajectory numbers
    Yfull=postfault.Yfull_mod;
    norm_min = 0;
    Normtt=zeros(n_itermax*n_MGPtrajmax,1);  % for observation: all Norm from each trajectory
     f1=evalin('base','f1');
     f2=evalin('base','f2');
%% Initialization
    flag_MGP=0;
    deltac_start=deltac_escape';
    theta_start=theta_escape';
    voltage_start=voltage_escape';
    plot(postfault.SEP_delta(2),postfault.SEP_delta(3),'ob','LineWidth',1.5,'MarkerSize',8); hold on;
%% MGP calculation
while(flag_MGP==0)
    [deltac_iter,theta_iter,voltage_iter,Normp,no_MGP,flag_MGP]=Fun_Cal_MGP_singletraj_SPM(deltac_start,theta_start,voltage_start,Tunit,n_itermax,norm_Tol,Yfull,preset);
    n_MGPtraj=n_MGPtraj+1;  % counter for iteration ++
    for i=1:size(Normp,1)
        noNorm=(n_MGPtraj-1)*n_itermax+i;
        Normtt(noNorm)=Normp(i);
    end
    if(flag_MGP==0)
         figure(f1);
         plot(deltac_iter(:,2),deltac_iter(:,3),'k-','LineWidth',1.5);
%         figure(f2);
%         plot((n_MGPtraj-1)*n_itermax+1:n_MGPtraj*n_itermax,Normp,'k-','LineWidth',1.5);
%         if(n_MGPtraj>1)
%             plot((n_MGPtraj-1)*n_itermax:(n_MGPtraj-1)*n_itermax+1,[norm_min Normp(1)],':r','LineWidth',1.5);
%         end
        norm_min = Normp(end);
    else
         figure(f1);
         plot(deltac_iter(1:no_MGP,2),deltac_iter(1:no_MGP,3));
        deltac_MGP=deltac_iter(no_MGP,:);
        num_Traj=n_MGPtraj;
        norm_min = Normp(no_MGP);
%         figure(f2);
%         plot((n_MGPtraj-1)*n_itermax+1:(n_MGPtraj-1)*n_itermax+no_MGP+1,Normp(1:no_MGP+1),'k-','LineWidth',1.5);
        break;
    end
  %  hold on;
    %% No MGP found in last iteartion process, update start point
    deltac_last=deltac_iter(n_itermax,:)';
    theta_last=theta_iter(n_itermax,:)';
    voltage_last=voltage_iter(n_itermax,:)';
    [deltac_update,theta_update,voltage_update,flag_update]=Fun_Cal_UpdateStartPoint_SPM(deltac_last,theta_last,voltage_last,preset,postfault);
    if(flag_update==1)
        deltac_starthis=deltac_start;
        deltac_start=deltac_update;
        theta_start=theta_update;
        voltage_start=voltage_update;
        if(norm(deltac_start-deltac_starthis)<1e-3)
            flag_MGP=1;
            deltac_MGP=deltac_update';
            num_Traj=n_MGPtraj;
            fprintf('MGP found since the iteration process reached a repeated status!\n');
            break;
        end
        if(norm(deltac_start-deltac_starthis)>0.5*norm(deltac_starthis-postfault.SEP_delta))
            deltac_start=deltac_last;
            theta_start=theta_last;
            voltage_start=voltage_last;
            fprintf('No update makes this time\n');
        end
        figure(f1);
        plot([deltac_last(2),deltac_update(2)],[deltac_last(3),deltac_update(3)],':k','LineWidth',1.5);
    else
        deltac_start=deltac_last;
        theta_start=theta_last;
        voltage_start=voltage_last;
    end
    if(n_MGPtraj>n_MGPtrajmax)
        error('No MGP found in %d times',n_MGPtraj);
    end
end