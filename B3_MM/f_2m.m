function dfdt = f_2m(x)


delta12 = x(1);
omega12 = x(2);
omegasum = x(3);

G1 = evalin('base','G1');
G2 = evalin('base','G2');
G12 = evalin('base','G12');
B12 = evalin('base','B12');
E1 = evalin('base','E1');
E2 = evalin('base','E2');
Pm1 = evalin('base','Pm1');
Pm2 = evalin('base','Pm2');
H1 = evalin('base','H1');
H2 = evalin('base','H2');
D1 = evalin('base','D1');
D2 = evalin('base','D2');

Pe1 = E1^2*(G12+G1) + E1*E2*G12*cos(delta12) + E1*E2*B12*sin(delta12);
Pe2 = E2^2*(G12+G2) + E1*E2*G12*cos(delta12) - E1*E2*B12*sin(delta12);

dfdt(1) = omega12;% delta12
dfdt(2) = Pm1/H1 - Pm2/H2 - Pe1/H1 + Pe2/H2 - omega12*(D1/H1+D2/H2)/2 - omegasum*(D1/H1-D2/H2)/2;%omega12
dfdt(3) = Pm1/H1 + Pm2/H2 - Pe1/H1 - Pe2/H2 - omega12*(D1/H1-D2/H2)/2 - omegasum*(D1/H1+D2/H2)/2;%omegasum


dfdt = dfdt.';

end