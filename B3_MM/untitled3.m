deltac = x_all(2000,1:2)';
m=preset.m;
delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
deltacc = [delta1c deltac(1) deltac(2)];
theta_net = x_all(2000,3:8)';
voltage_net = x_all(2000,9:14)';
temp_ini = [theta_net;voltage_net];
system = "postfault";
ccc=Fun_AEfslove_SPM(temp_ini,deltacc,preset,system)



%% Jacobi matrix calculation
V_net=voltage_net;
delta_net=theta_net;
deltac=deltacc;


    G=real(postfault.Yfull_mod);    
    B=imag(postfault.Yfull_mod);
    J11=zeros((nbus-ngen),(nbus-ngen)); % Pnet/delta
    J12=zeros((nbus-ngen),(nbus-ngen)); % Pnet/V
    J21=zeros((nbus-ngen),(nbus-ngen)); % Qnet/delta
    J22=zeros((nbus-ngen),(nbus-ngen)); % Qnet/V
    for i=1:(nbus-ngen)
        J12(i,i)=2*V_net(i)*G(i+ngen,i+ngen);
        J22(i,i)=-2*V_net(i)*B(i+ngen,i+ngen);
        for j=1:ngen
            ddelta=delta_net(i)-deltac(j);
            J11(i,i)=J11(i,i)+V_net(i)*E(j)*B(i+ngen,j)*cos(ddelta)-V_net(i)*E(j)*G(i+ngen,j)*sin(ddelta);
            J12(i,i)=J12(i,i)+E(j)*B(i+ngen,j)*sin(ddelta)+E(j)*G(i+ngen,j)*cos(ddelta);
            J21(i,i)=J21(i,i)+V_net(i)*E(j)*B(i+ngen,j)*sin(ddelta)+V_net(i)*E(j)*G(i+ngen,j)*cos(ddelta);
            J22(i,i)=J22(i,i)-E(j)*B(i+ngen,j)*cos(ddelta)+E(j)*G(i+ngen,j)*sin(ddelta);
        end
        for j=1:(nbus-ngen)
            if (i~=j)
                ddelta=delta_net(i)-delta_net(j);
                J11(i,i)=J11(i,i)+V_net(i)*V_net(j)*B(i+ngen,j+ngen)*cos(ddelta)-V_net(i)*V_net(j)*G(i+ngen,j+ngen)*sin(ddelta);
                J11(i,j)=-V_net(i)*V_net(j)*B(i+ngen,j+ngen)*cos(ddelta)+V_net(i)*V_net(j)*G(i+ngen,j+ngen)*sin(ddelta);

                J12(i,i)=J12(i,i)+V_net(j)*B(i+ngen,j+ngen)*sin(ddelta)+V_net(j)*G(i+ngen,j+ngen)*cos(ddelta);
                J12(i,j)=V_net(i)*B(i+ngen,j+ngen)*sin(ddelta)+V_net(i)*G(i+ngen,j+ngen)*cos(ddelta);

                J21(i,i)=J21(i,i)+V_net(i)*V_net(j)*B(i+ngen,j+ngen)*sin(ddelta)+V_net(i)*V_net(j)*G(i+ngen,j+ngen)*cos(ddelta);
                J21(i,j)=-V_net(i)*V_net(j)*B(i+ngen,j+ngen)*sin(ddelta)-V_net(i)*V_net(j)*G(i+ngen,j+ngen)*cos(ddelta);

                J22(i,i)=J22(i,i)-V_net(j)*B(i+ngen,j+ngen)*cos(ddelta)+V_net(j)*G(i+ngen,j+ngen)*sin(ddelta);
                J22(i,j)=-V_net(i)*B(i+ngen,j+ngen)*cos(ddelta)+V_net(i)*G(i+ngen,j+ngen)*sin(ddelta);
            end
        end
    end
    J=[J11, J12;J21, J22]
    eig(J)