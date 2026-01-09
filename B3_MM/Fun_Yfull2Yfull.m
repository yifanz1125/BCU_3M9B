%% function: modify structure preserved admittance matrix
function [Y_full_modified,Transform]=Fun_Yfull2Yfull(Y_full,pfdata,faultflag)
%% block matrix
    if(faultflag==0)
        nbus=pfdata.bus.num;
        no_faultbus=nbus+1;
    else
        nbus=pfdata.bus.num-1;
        no_faultbus=faultflag(2);
    end
    ngen=pfdata.bus.numgen;
    Ynn=zeros(ngen,ngen);
    Ynr=zeros(ngen,nbus-ngen);
    Yrn=zeros(nbus-ngen,ngen);
    Yrr=zeros(nbus-ngen,nbus-ngen);    
    % Ynn
        for i=1:ngen      
            for j=1:ngen
                no_genself=pfdata.gen.no(i,1);
                no_geninte=pfdata.gen.no(j,1);
                if(no_genself>no_faultbus)
                    no_genself=no_genself-1;
                end
                if(no_geninte>no_faultbus)
                    no_geninte=no_geninte-1;
                end
                Ynn(i,j)=Y_full(no_genself,no_geninte);
            end
        end
    % Ynr
        no_withoutgen=zeros(nbus-ngen,1);
        k=1;
        for i=1:nbus
            flag_gen=0;
            for j=1:ngen
                genno = pfdata.gen.no(j);
                if(genno>no_faultbus)
                    genno = genno-1;
                end
                if(i==genno)
                    flag_gen=1;
                end
            end
            if(flag_gen==0)
                no_withoutgen(k,1)=i;
                k=k+1;
            end
            flag_gen=0; 
        end
        for i=1:ngen      
            for j=1:(nbus-ngen)
                no_genself=pfdata.gen.no(i,1);
                if(no_genself>no_faultbus)
                    no_genself=no_genself-1;
                end
                no_loadinte=no_withoutgen(j,1);
                Ynr(i,j)=Y_full(no_genself,no_loadinte);
            end
        end    
    % Yrn
        for i=1:(nbus-ngen)      
        for j=1:ngen
            no_loadself=no_withoutgen(i,1);
            no_geninte=pfdata.gen.no(j,1);
            if(no_geninte>no_faultbus)
                no_geninte=no_geninte-1;
            end
            Yrn(i,j)=Y_full(no_loadself,no_geninte);
        end
        end   
    % Yrr
        for i=1:(nbus-ngen)      
        for j=1:(nbus-ngen) 
            no_loadself=no_withoutgen(i,1);
            no_loadinte=no_withoutgen(j,1);
            Yrr(i,j)=Y_full(no_loadself,no_loadinte);
        end
        end    
%% modified matrix
    Y_full_modified=[Ynn,Ynr;Yrn,Yrr];
    k=1;
     for i=1:pfdata.bus.num
            flag_gen=0;
            for j=1:ngen
                genno = pfdata.gen.no(j);
                if(i==genno)
                    flag_gen=1;
                end
            end
            if(flag_gen==0 && i~=no_faultbus)
                no_withoutgen(k,1)=i;
                k=k+1;
            end
    end
    Transform = [pfdata.gen.no;no_withoutgen];
end
