%% Calculate the internal voltage E∠δ of generators, based on the power flow results
function EMF=Fun_Cal_GenEMF(flagxd,pfdata,xd1)
    Sbase=pfdata.Sbase;
%% E-delta calculation
    EMF=zeros(pfdata.bus.numgen,2);
    if(flagxd==1)
        for i=1:pfdata.bus.numgen
            Vpu=pfdata.gen.voltage(i,1);
            Ppu=pfdata.gen.PQ(i,1)/Sbase;
            Qpu=pfdata.gen.PQ(i,2)/Sbase;
            EMF(i,1)=sqrt((Vpu+Qpu*xd1(i)/Vpu)^2+(Ppu*xd1(i)/Vpu)^2);
            EMF(i,2)=asin((Ppu*xd1(i)/Vpu)/(Vpu+Qpu*xd1(i)/Vpu));
        end
        EMF(:,2)=EMF(:,2)+pfdata.gen.voltage(:,2);  % unit: rad
    else
        EMF(:,1)=pfdata.gen.voltage(:,1);
        EMF(:,2)=pfdata.gen.voltage(:,2);  % unit: rad
    end    
end