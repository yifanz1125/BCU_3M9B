%% unit: omega in rad/s (in SRF)
function f= Fun_SEPfslove(delta_omega,preset,state,basevalue)
    Y_red=state.Yred;
    G=real(Y_red);    
    B=imag(Y_red);
    ngen=size(Y_red,1);
    m=preset.m;
    d=preset.d;
    Pm=preset.Pmpu;
    E=preset.Epu;
    omegab=basevalue.omegab;
    

% f=State.Pm*0;
    delta=zeros(ngen,1);
    Pe=zeros(ngen,1);
        omega=delta_omega(ngen);
        for i=1:ngen-1
            delta(i)=delta_omega(i);
        end

    % Pe calculation
        for i=1:ngen
            for j=1:ngen
                ddelta=delta(i,1)-delta(j,1);
                Pe(i)=Pe(i)+E(i)*E(j)*B(i,j)*sin(ddelta)+E(i)*E(j)*G(i,j)*cos(ddelta);
            end
        end

        PCOI=sum(Pm-Pe)-sum(d)*omega;

    for i=1:ngen-1
        f(i)=Pm(i)-Pe(i)-m(i)/sum(m)*PCOI-d(i)*omega;
    end
    f(ngen)= sum(Pm-Pe)-sum(d)*omega;    % the nth equ

end

