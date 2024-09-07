

% system = "postfault";
% delta0=x_fault_all(end,1:3)';
% omega0=x_fault_all(end,16:18)';
% %delta_net0=x_fault_all(end,4:8)';
% %delta_net0=[delta_net0(1:(fault.faultbus-ngen-1)); 0; delta_net0((fault.faultbus-ngen):end)];
% %voltage_net0=x_fault_all1(end,10:14)';
% %voltage_net0=[voltage_net0(1:(fault.faultbus-ngen-1)); 1; voltage_net0((fault.faultbus-ngen):end)];
% %[delta_net_s,Voltage_net_s,flag_iter,n_iter,err] = Fun_AEfslove(zeros(6,1),ones(6,1),delta0,preset,Basevalue,1e4,1e-10);
% %[t_postfault, x_postfault_all] = ode15s(@f_timedomain_DAE_withflove,[Iter.Trecover,Iter.Ttotal],[delta0; delta_net_s; Voltage_net_s; omega0],options);
% [t_postfault, x_postfault_all] = ode15s(@f_timedomain_DAE_withflove,[Iter.Trecover,Iter.Ttotal],[delta0; omega0],odeset('RelTol',1e-10,'AbsTol',1e-8*ones(1,6)));
% 
% function dfdt = f_timedomain_DAE(t,x)
%     dfdt = F_3M9B_SP_DAE(x);
% end
% function dfdt = f_timedomain_DAE_withflove(t,deltaomega)
%     preset = evalin('base','preset');
%     Results_fsolve=fsolve(@(x)Fun_AEfslove_SPM(x,deltaomega(1:3),preset),[zeros(6,1);ones(6,1)],optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-9));
%     xx = [deltaomega(1:3);Results_fsolve;deltaomega(4:6)];
%     dfdtt = F_3M9B_SP_DAE(xx);
%     dfdt=[dfdtt(1:3);dfdtt(16:18)];
% end




system = "fault1";
[delta_net_s,V_net_s,flag_iter,n_iter,err] = Fun_AEiteration_SPM(delta_net0,voltage_net0,delta0,preset,Basevalue,1e5,1e-10)
Results_fsolve=fsolve(@(x)Fun_AEfslove_SPM(x,delta0,preset),[delta_net_s(1:5);V_net_s(1:5)],optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','iter','TolX',1e-9))

Fun_AEfslove_SPM([delta_net_s(1:5);V_net_s(1:5)],delta0,preset)