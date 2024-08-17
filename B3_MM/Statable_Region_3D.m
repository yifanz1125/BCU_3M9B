clear;

%% Parameter
Z1 = 0.05+0.5j;
Z2 = 0.01+0.3j;
Zl = 0.2844+0.0306j;
Pm1 = 1.33;
Pm2 = 0.6;
H1 = 0.5;
H2 = 0.5;
D1 = 0.4;
D2 = 0.5;
Y12 = 1/(Z1+Z2+Z1*Z2/Zl);
Y1 = 1/(Z1+Zl+Z1*Zl/Z2);
Y2 = 1/(Z2+Zl+Z2*Zl/Z1);
G1 = real(Y1);
G2 = real(Y2);
G12 = -real(Y12);
B12 = -imag(Y12);
E1 = 1; E2 = 1;

%%
global N;
try
    N;
catch
    N = 2;
end
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(N,n);
for mm = 1:n
   x_set(1,mm) = x(mm);
end
torralence = 1e-2; 
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f1,xep,options);

    %xep = mod(xep,2*pi);
    
    if maxabs(ferr) < torralence
        if isnewxep(ep_set,xep,torralence)
           
            [V,Lambda]=eig(A);
            Lambda = diag(Lambda);
            sig = sign(sign(real(Lambda))+0.1); % zero counted as positive
            sig = (sig + 1)/2;                  % [0,1]
            flag = sum(sig);                    % number of non-negative eigenvalues

            v = V(:,~sig);                      % the stable sub-space
            
            ep_set(m).xep = real(xep); %#ok<*SAGROW> 
            ep_set(m).A = A;
            ep_set(m).Lambda = Lambda;
            ep_set(m).V = V;   
            ep_set(m).v = v;     % stable eigenvectors of unstable ep 
            ep_set(m).flag = flag;
           
            m = m+1;
            if flag == 0
            jacob=A;
            xeps=xep;
            end
        end
    end
end

for m = 1:length(ep_set)
    disp_v('Index',m);
    disp_v('Equilibrium',ep_set(m).xep);
    disp_v('Eigenvalue', ep_set(m).Lambda);
    disp_v('Eigenvector',ep_set(m).V);
end
if N == 2
for n = 1:length(ep_set)
        m = (n-1)*3;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
    
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
    
        ep_set_ext(m+3).xep(1) = ep_set(n).xep(1) + 2*pi;
    
end
figure;
hold on;
grid on;
color_code = {'blue','magenta','red','black'};
axis([-2*pi,2*pi,-2*pi,2*pi]);
for m = 1 : length(ep_set_ext)
        xep = ep_set_ext(m).xep;
        flag= ep_set_ext(m).flag;
        scatter(xep(1),xep(2),color_code{flag+1});
       
        if flag == 1
            v = ep_set_ext(m).v;
            perturb = 1e-2;
            [~ , x_p] = ode45(@f_backward,[0,50],xep+v*perturb,odeset('RelTol',1e-5));
            [~ , x_n] = ode45(@f_backward,[0,50],xep-v*perturb,odeset('RelTol',1e-5));
            x_all = [flip(x_n,1);x_p];
            plot(x_all(:,1),x_all(:,2),'k-','linewidth',1.5);
        end        
end
elseif N == 3
for n = 1:length(ep_set)
        m = (n-1)*3;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+3).xep(1) = ep_set(n).xep(1) + 2*pi;   
end
figure;

color_code = {'blue','magenta','green','red','cyan','blue','magenta','green','red','cyan'};
for m = 1:length(ep_set_ext)
    xep = ep_set_ext(m).xep;
    flag = ep_set_ext(m).flag;
    scatter3(xep(1),xep(2),xep(3),color_code{flag+1});
    xlabel('delta12');
    ylabel('omega12');
    zlabel('omega_{sum}');
    hold on;
end    
n = 1;
for m = 1 : length(ep_set_ext)        
    flag = ep_set_ext(m).flag;
    if flag 
        xep = ep_set_ext(m).xep;
        v = ep_set_ext(m).v;
        perturb = 5e-1;

        if flag == 1
            for alpha = (0:0.005:1)*2*pi
                vp = v(:,1)*sin(alpha) + v(:,2)*cos(alpha);
                [~ , x_all] = ode45(@f_backward,[0,10],xep+vp*perturb);
                plot3(x_all(:,1),x_all(:,2),x_all(:,3),color_code{n});
            end
            n = n + 1;
        elseif flag == 2
            for beta = [-1,1]
                vp = v*beta;
                [~ , x_all] = ode45(@f_backward,[0,2],xep+vp*perturb);
                plot3(x_all(:,1),x_all(:,2),x_all(:,3),color_code{n});
            end
            n = n + 1;
        end
        
    end
end
% plot trajectory
[tt , x_all] = ode78(@f_forward,[0,200],[-1.637, -6.87, -20],odeset('RelTol',1e-5));
[tt2 , x_all2] = ode78(@f_forward,[0,200],[0.311-4*pi, 4.75, 10],odeset('RelTol',1e-5));
[tt3 , x_all3] = ode78(@f_forward,[0,200],[-1.16, -9.67, -20],odeset('RelTol',1e-5));
[tt4 , x_all4] = ode78(@f_forward,[0,200],[1.02, 2.11, 10],odeset('RelTol',1e-5));
[tt5 , x_all5] = ode78(@f_forward,[0,200],[1.02, 1.11, -20],odeset('RelTol',1e-5));
[tt6 , x_all6] = ode78(@f_forward,[0,1200],[-2.667, 0.045+1, -60],odeset('RelTol',1e-5));%[-2.667, 0.045+1, 80]
[tt7 , x_all7] = ode78(@f_forward,[0,200],[-1, -1.807, 20],odeset('RelTol',1e-5));
[tt8 , x_all8] = ode78(@f_forward,[0,200],[3, 0, 0],odeset('RelTol',1e-5));
plot3(x_all(:,1),x_all(:,2),x_all(:,3),'black','linewidth',1.5);
plot3(x_all2(:,1),x_all2(:,2),x_all2(:,3),'black','linewidth',1.5);
plot3(x_all3(:,1),x_all3(:,2),x_all3(:,3),'black','linewidth',1.5);
plot3(x_all4(:,1),x_all4(:,2),x_all4(:,3),'black','linewidth',1.5);
plot3(x_all5(:,1),x_all5(:,2),x_all5(:,3),'black','linewidth',1.5);
plot3(x_all6(:,1),x_all6(:,2),x_all6(:,3),'black','linewidth',1.5);
plot3(x_all7(:,1),x_all7(:,2),x_all7(:,3),'black','linewidth',1.5);
plot3(x_all8(:,1),x_all8(:,2),x_all8(:,3),'black','linewidth',1.5);
axis([-2*pi 2*pi -15 10 -100 100]);
end





function yes = isnewxep(ep_set,xep,torr)
    if isempty(ep_set)
        yes = 1;
        return;
    end
    minerr = inf;
    for m = 1 : length(ep_set)
        err = abs(xep - ep_set(m).xep);
        err = min(err, abs(2*pi-err));
        err = max(err);
        if minerr > err
            minerr = err;
        end
    end
    if(minerr>torr)
        yes = 1;
    else
        yes = 0;
    end
end
function disp_v(msg,v)
    disp([msg '=']);
    disp(v);
end
function dfdt = f1(x)  
 global N; 
    switch N
        case 2
            dfdt = f_2m_reduce(x);
        case 3
            dfdt = f_2m(x);
    end 
end
function dfdt = f_forward(t,x)
    dfdt = f1(x);
end

function dfdt = f_backward(t,x)
    dfdt = -f1(x);
end


