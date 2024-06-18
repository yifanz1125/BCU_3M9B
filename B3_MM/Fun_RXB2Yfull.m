%% function: transfer RXB (from matpower results) to admittance matrix (full order)
function Yfull=Fun_RXB2Yfull(RXB,pfdata)
    Yii0=zeros(pfdata.bus.num,2);
    Yij_re=zeros(pfdata.bus.num,pfdata.bus.num);
    Yij_im=zeros(pfdata.bus.num,pfdata.bus.num);
    for i=1:pfdata.bus.num  % search for the bus i
        for k=1:size(RXB,1)  % scan the k line
            if(RXB(k,1)==i||RXB(k,2)==i)
                Yii0(i,1)=Yii0(i,1)+RXB(k,3)/(RXB(k,3)^2+RXB(k,4)^2);
                Yii0(i,2)=Yii0(i,2)+RXB(k,5)/2-RXB(k,4)/(RXB(k,3)^2+RXB(k,4)^2);
                if(RXB(k,1)==i)
                    j=RXB(k,2);
                else
                    j=RXB(k,1);
                end
                Yij_re(i,j)=-1*RXB(k,3)/(RXB(k,3)^2+RXB(k,4)^2);
                Yij_im(i,j)=RXB(k,4)/(RXB(k,3)^2+RXB(k,4)^2);
                Yij_re(j,i)=Yij_re(i,j);
                Yij_im(j,i)=Yij_im(i,j);
            end
        end
    end
    Yij=complex(Yij_re,Yij_im);        
    Yii=complex(Yii0(:,1),Yii0(:,2));
    Yii_mat=diag(Yii);
    Yfull=Yij+Yii_mat;
end