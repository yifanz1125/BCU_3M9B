m=1;
figure;
hold on;
grid on;
color_code = {'blue','magenta','red','black'};
axis([-2*pi,2*pi,-2*pi,2*pi]);
xep = ep_set_ext(m).xep;
flag= ep_set_ext(m).flag_reduce;
scatter(xep(1),xep(2),color_code{flag+1});
       
v = ep_set_ext(m).v_reduce;
perturb = 1e-2;
perturb_xep_p = xep+v*perturb;
perturb_xep_n = xep-v*perturb;
theta_net = ep_set_ext(m).delta_net_ep;
voltage_net = ep_set_ext(m).voltage_net_ep;
temp_ini = [theta_net;voltage_net];
delta1c= -mmm(2:ngen)'*perturb_xep_p/mmm(1);
deltacc_p = [delta1c perturb_xep_p(1) perturb_xep_p(2)];
delta1c= -mmm(2:ngen)'*perturb_xep_n/mmm(1);
deltacc_n = [delta1c perturb_xep_n(1) perturb_xep_n(2)];
[results_p,fval,exitflag,output,jacobian]=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc_p,preset,"postfault"),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
perturb_p = [perturb_xep_p;results_p];
[results_n,fval,exitflag,output,jacobian]=fsolve(@(x)Fun_AEfslove_SPM(x,deltacc_n,preset,"postfault"),temp_ini,optimset('TolFun',1e-20,'MaxFunEvals',1e5,'Maxiter',1e5,'Display','off','TolX',1e-5));
perturb_n = [perturb_xep_n;results_n];

Fun_AEfslove_SPM(results_p,deltacc_p,preset,"postfault")
f_reducedstate_SPM(perturb_p)
Fun_AEfslove_SPM(results_n,deltacc_n,preset,"postfault")
f_reducedstate_SPM(perturb_n)

%%
M = diag([ones(2,1); ones(12,1)*1e-15]);
options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,2),1e-12*ones(1,12)]);
[~ , x_p] = ode15s(@f_backward,[0,10],perturb_p,options);
[~ , x_n] = ode15s(@f_backward,[0,10],perturb_n,options);
x_all = [flip(x_n,1);x_p];
plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);


function dfdt = f_backward(t,x)
    dfdt = f_reducedstate_SPM_backward(x);
end