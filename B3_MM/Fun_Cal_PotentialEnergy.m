%% function: calculate potential energy (defined in original system) 
%% Input: thetac_end(n*1)--current point thetac_start--used for path estimation
%% Output: Ep1--equivalent Pm Ep2--equivalent Bij*sinδ item Ep2--path-relevant (lossy network)
function [Ep1,Ep2,Ep3]=Fun_Cal_PotentialEnergy(preset,postfault,thetac_start,thetac_end)
    thetac_SEP=postfault.SEP_delta;
    Y_red=postfault.Yred;
    ngen=size(Y_red,1);
    G_post=real(Y_red);
    B_post=imag(Y_red);
    Pm=preset.Pmpu;
    E=preset.Epu;
    i=preset.m;
    d=preset.d;
    Ep1=0;
    Ep2=0;
    Ep3=0;
    Ep3_rad=0;
    for i=1:ngen
        Ep1=Ep1+(-1)*(Pm(i)-E(i)^2*G_post(i,i))*(thetac_end(i,1)-thetac_SEP(i,1));
    end
    for i=1:ngen-1
        for j=i+1:ngen
            Ep2=Ep2+(-1)*E(i)*E(j)*B_post(i,j)*(cos(thetac_end(i,1)-thetac_end(j,1))-cos(thetac_SEP(i,1)-thetac_SEP(j,1)));
        end
    end
% Path-relevant Potential Energy Estimation
% radial one-as ref
        for i=1:ngen-1 % 耗散势能
            for j=i+1:ngen
                dtheta_i=thetac_end(i) -thetac_start(i);
                dtheta_j=thetac_end(j) -thetac_start(j);
                dtheta_ij=dtheta_i-dtheta_j;
                if(abs(dtheta_ij)>1e-6) % 避免数值病态问题
                    adb=(dtheta_i+dtheta_j)/dtheta_ij;                
                else
                    adb=dtheta_i+dtheta_j;
                end
                Ep3_rad= Ep3_rad+E(i)*E(j)*G_post(i,j)*adb*(sin(thetac_end(i)-thetac_end(j))-sin(thetac_start(i)-thetac_start(j)));
            end
        end



    if(preset.PathEnergyCal==0)
        for i=1:ngen-1 % 耗散势能
            for j=i+1:ngen
                dtheta_i=thetac_end(i) -thetac_start(i);
                dtheta_j=thetac_end(j) -thetac_start(j);
                dtheta_ij=dtheta_i-dtheta_j;
                if(abs(dtheta_ij)>1e-7) % 避免数值病态问题
                    adb=(dtheta_i+dtheta_j)/dtheta_ij;                
                else
                    adb=dtheta_i+dtheta_j;
                end
                Ep3= Ep3+E(i)*E(j)*G_post(i,j)*adb*(sin(thetac_end(i)-thetac_end(j))-sin(thetac_start(i)-thetac_start(j)));
            end
        end
    elseif(preset.PathEnergyCal==-1)
        Ep3=0;
    else
        n_mid=preset.PathEnergyCal-1;   % numbers of inserted mid-point (n_mid=0--trapezoidal from start to end)  
        dtheta=thetac_end-thetac_start;  % dtheta(i)=theta_end(i)-theta_start(i)
        unit_dtheta=dtheta/(n_mid+1);
        for i=1:ngen
            for j=1:ngen
                if(i~=j)
                    for m=1:n_mid+1
                    Ep3=Ep3+E(i)*E(j)*G_post(i,j)*0.5*unit_dtheta(i)...,
                        *(cos(thetac_start(i)+(m-1)*unit_dtheta(i)-thetac_start(j)-(m-1)*unit_dtheta(j))+...,
                        cos(thetac_start(i)+m*unit_dtheta(i)-thetac_start(j)-m*unit_dtheta(j)));
                    end
                end
            end
        end
    end

end