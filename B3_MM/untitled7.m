G_post=real(Yfull);
    B_post=imag(Yfull);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    ngen=preset.ngen;
    nbus=preset.nbus;
    %% Initialization
    flag_normmin=0;
    Pe=zeros(n_itermax,ngen);
    Pcoi=zeros(n_itermax,1);
    deltac_iter=zeros(n_itermax,ngen);
    theta_iter=zeros(n_itermax,nbus-ngen);
    voltage_iter=zeros(n_itermax,nbus-ngen);
    Normp=zeros(n_itermax,1);
    norm_min=0;
    no_MGP=0;
    

    %% iteration process
    M = diag([ones(2,1); ones(12,1)*1e-15]);
    options = odeset('Mass',M,'RelTol',1e-10,'AbsTol',[1e-8*ones(1,2),1e-12*ones(1,12)]);
    [theta_start,voltage_start,flag_iter,n_iter,err] = Fun_AEiteration_SPM(theta_start,voltage_start,deltac_start,preset,Basevalue,system,1e4,1e-10);
    [t_postfault, x_postfault_all] = ode15s(@fred,0:Tunit:Tunit*(n_itermax-1),[deltac_start(2:3); theta_start; voltage_start],options);
    deltac_temp = x_postfault_all(:,1:2);
    delta1c= -m(2:ngen)'*deltac_temp'/m(1);
    deltac_iter = [delta1c' deltac_temp];
    theta_iter = x_postfault_all(:,3:8);
    voltage_iter = x_postfault_all(:,9:14);



    %% 
    for tm=1:n_itermax
        % Pe calculation
        for  i=1:ngen
            for j=1:ngen
                ddelta=deltac_iter(tm,i)-deltac_iter(tm,j);
                Pe(tm,i)=Pe(tm,i)+E(i)*E(j)*(G_post(i,j)*cos(ddelta)+B_post(i,j)*sin(ddelta));
            end
            for l=1:(nbus-ngen)
                ddelta=deltac_iter(tm,i)-theta_iter(tm,l);
                Pe(tm,i)=Pe(tm,i)+E(i)*voltage_iter(tm,l)*B_post(i,l+ngen)*sin(ddelta)+E(i)*voltage_iter(tm,l)*G_post(i,l+ngen)*cos(ddelta);
            end
        end
        Pcoi(tm)=sum(Pm)-sum(Pe(tm,:));

        
        Normp(tm)=norm((Pm'-Pe(tm,:)-Pcoi(tm,1)/sum(m)*m'));%./d'
        if(tm==1)
            norm_min=Normp(1);
        else
            if(Normp(tm)<norm_min)
                norm_min=Normp(tm);
            end
            if(Normp(tm)-Normp(tm-1)>norm_Tol&&Normp(tm-1)==norm_min&&tm~=2&&norm_min<1e-1)
                flag_normmin=1;
                no_MGP=tm-1;
                fprintf('MGP found since the norm is rather small and a minimum found');
                break;
            end
        end
    end
%%

%%
function dfdt = fred(t,x)
    dfdt = f_reducedstate_SPM(x);
end