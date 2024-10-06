function dfdt = f_reducedstate(deltac)

postfault = evalin('base','postfault');
preset = evalin('base','preset');

delta2c=deltac(1);
delta3c=deltac(2);

Y_red=postfault.Yred;
G=real(Y_red);    
B=imag(Y_red);
ngen=size(Y_red,1);
m=preset.m;
d=preset.d;
Pm=preset.Pmpu;
E=preset.Epu;
mT=sum(m,1);
Pe=zeros(1,ngen);

delta1c= -m(2:ngen)'*deltac(1:2)/m(1);
deltacc = [delta1c delta2c delta3c];

% Pe calculation
for i=1:ngen
    for j=1:ngen
        ddelta=deltacc(i)-deltacc(j);
        Pe(i)=Pe(i)+E(i)*E(j)*(G(i,j)*cos(ddelta)+B(i,j)*sin(ddelta));
    end
end
Pcoi=sum(Pm)-sum(Pe);

dfdt(1) =  (Pm(2)-Pe(2)-Pcoi/sum(m)*m(2)); 
dfdt(2) =  (Pm(3)-Pe(3)-Pcoi/sum(m)*m(3));


dfdt = dfdt.';

end