function  [thetac_iter,Normp,no_MGP,flag_normmin]=Fun_Cal_MGP_singletraj(thetac_start,Tunit,n_itermax,norm_Tol,Yred_post,preset)
    %% Settings
    G_post=real(Yred_post);
    B_post=imag(Yred_post);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    ngen=size(Yred_post,1);
    %% Initialization
    flag_normmin=0;
    Pe=zeros(n_itermax,ngen);
    Pcoi=zeros(n_itermax,1);
    thetac_iter=zeros(n_itermax,ngen);
    Normp=zeros(n_itermax,1);
    norm_min=0;
    no_MGP=0;
    thetac_iter(1,:)=thetac_start;
    thetac_tmp=zeros(1,ngen);
    %% iteration process
    for tm=1:n_itermax
        % Pe calculation
        for i=1:ngen
            for j=1:ngen
                ddelta=thetac_iter(tm,i)-thetac_iter(tm,j);
                Pe(tm,i)=Pe(tm,i)+E(i)*E(j)*B_post(i,j)*sin(ddelta)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
            end
        end
        Pcoi(tm)=sum(Pm)-sum(Pe(tm,:));
        if(tm<n_itermax)
            for i=1:ngen-1
                %thetac_tmp(i)=thetac_iter(tm,i)-thetac_iter(tm,ngen)+((Pm(i)-Pe(tm,i)-Pcoi(tm,1)/sum(m)*m(i))/d(i))*Tunit;
                thetac_tmp(i)=thetac_iter(tm,i)-thetac_iter(tm,ngen)+((Pm(i)-Pe(tm,i)-Pcoi(tm,1)/sum(m)*m(i))/d(i)-(Pm(ngen)-Pe(tm,ngen)-Pcoi(tm,1)/sum(m)*m(ngen))/d(ngen))*Tunit;
            end
            thetac_tmp(ngen)=0;
            thetac_iter(tm+1,:)=thetac_tmp-m'*thetac_tmp'/sum(m);
        end

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
end
