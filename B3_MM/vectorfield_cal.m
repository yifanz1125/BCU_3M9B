%% unit: omega in rad/s (in SRF)
function f= vectorfield_cal(thetac,Yred_post,preset)
    %% Settings
    G_post=real(Yred_post);
    B_post=imag(Yred_post);
    Pm=preset.Pmpu;
    E=preset.Epu;
    m=preset.m;
    d=preset.d;
    ngen=size(Yred_post,1);
    Pe=zeros(1,ngen);

    thetac= thetac;
    for i=1:ngen
            for j=1:ngen
                ddelta=thetac(i)-thetac(j);
                Pe(i)=Pe(i)+E(i)*E(j)*B_post(i,j)*sin(ddelta)+E(i)*E(j)*G_post(i,j)*cos(ddelta);
            end
    end
    Pcoi=sum(Pm)-sum(Pe);
    Normp=norm((Pm'-Pe-Pcoi/sum(m)*m'));
    f=Normp;

end

