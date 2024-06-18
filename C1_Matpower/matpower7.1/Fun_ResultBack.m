% output the powerflow result back to main function
function pfdata=Fun_ResultBack(casename)
    folder_matpower='C:\Users\yz7521\OneDrive - Imperial College London\BCU\C1_Matpower\matpower7.1';
    folder_matpower2='C:\Users\yz7521\OneDrive - Imperial College London\BCU\C1_Matpower\matpower7.1\lib';
    folder_matpower3='C:\Users\yz7521\OneDrive - Imperial College London\BCU\C1_Matpower\matpower7.1\data';
    addpath(genpath(folder_matpower));  addpath(genpath(folder_matpower2));  addpath(genpath(folder_matpower3));
    install_matpower(1,0,0,1);
    results=runpf(casename);
    
    if(results.success)
        pfdata.voltage=results.bus(:,8:9);
        pfdata.bus.numgen=size(results.gen,1);  % number of generators
        pfdata.gen.no=results.gen(:,1); % generators bus no
        pfdata.gen.PQ=results.gen(:,2:3);   %generators PQ output
        pfdata.gen.voltage=zeros(pfdata.bus.numgen,2);  % generator bus voltages (mag and angle)
        for i=1:pfdata.bus.numgen
            genno=pfdata.gen.no(i,1);
            pfdata.gen.voltage(i,1:2)=pfdata.voltage(genno,1:2);
        end
    
        pfdata.Sbase=results.baseMVA;
        pfdata.bus.PQ=results.bus(:,3:4);
        pfdata.bus.num=size(results.bus,1); % number of network buses
        pfdata.bus.numload=0;   % number of PQ output buses
        loadno=zeros(pfdata.bus.num,1);
        loadPQ=zeros(pfdata.bus.num,2);    
        for i=1:pfdata.bus.num
            Pbus=results.bus(i,3);
            Qbus=results.bus(i,4);
            if(Pbus~=0||Qbus~=0)
                pfdata.bus.numload=pfdata.bus.numload+1;
                loadno(pfdata.bus.numload,1)=results.bus(i,1);
                loadPQ(pfdata.bus.numload,1)=Pbus;
                loadPQ(pfdata.bus.numload,2)=Qbus;           
            end
        end
    
        pfdata.load.no=zeros(pfdata.bus.numload,1); pfdata.load.PQ=zeros(pfdata.bus.numload,2);
        pfdata.load.no=loadno(1:pfdata.bus.numload,1);  % PQ output bus no
        pfdata.load.PQ=loadPQ(1:pfdata.bus.numload,:);  % PQ output bus PQ
        for i=1:pfdata.bus.numload
            loadno=pfdata.load.no(i,1);
            pfdata.load.voltage(i,1:2)=pfdata.voltage(loadno,1:2);
        end
        pfdata.gen.voltage(:,2)=pfdata.gen.voltage(:,2)*pi/180;
        pfdata.load.voltage(:,2)=pfdata.load.voltage(:,2)*pi/180;
        pfdata.branch.RXB=results.branch(:,1:5);% fbus	tbus	r	x	b
        pfdata.branch.powerflow=zeros(size(results.branch,1),6);
        pfdata.branch.powerflow(:,1:2)=results.branch(:,1:2);   % FromBus no          TobusInjection no
        pfdata.branch.powerflow(:,3:6)=results.branch(:,14:17); % FromBusInjection PQ TobusInjection PQ
    else
        fprinf('Powerflow calculation failed!');
    end
end