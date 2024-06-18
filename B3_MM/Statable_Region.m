N=2;
x=(0:0.1:1)*2*pi;
n = length(x);
x_set = zeros(N,n^N);
for m = 0:(n^N - 1)
    for l = 1:N
        index = floor(m/n^(l-1));
        index = mod(index,n)+1;
        x_set(l,m+1) = x(index);
    end
end
torralence = 1e-2; 
%% calculate EPs
m = 1;
ep_set = [];
options = optimoptions('fsolve','FunctionTolerance',1e-10,'MaxIterations',10000,'OptimalityTolerance',1e-10);
for n = 1:length(x_set(1,:))
    xep = x_set(:,n);
    [xep,ferr,~,~,A] = fsolve(@f,xep,options);

    xep = mod(xep,2*pi);
    
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


for n = 1:length(ep_set)
        m = (n-1)*6;
        ep_set_ext(m+1)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+2)=ep_set(n);
        ep_set_ext(m+3)=ep_set(n);
        ep_set_ext(m+4)=ep_set(n);
        ep_set_ext(m+5)=ep_set(n); %#ok<*AGROW> 
        ep_set_ext(m+6)=ep_set(n);
    
        ep_set_ext(m+2).xep(1) = ep_set(n).xep(1) - 2*pi;
    
        ep_set_ext(m+3).xep(2) = ep_set(n).xep(2) - 2*pi;
    
        ep_set_ext(m+4).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+4).xep(2) = ep_set(n).xep(2) - 2*pi;
        ep_set_ext(m+5).xep(1) = ep_set(n).xep(1);
        ep_set_ext(m+5).xep(2) = ep_set(n).xep(2) - 4*pi;

        ep_set_ext(m+6).xep(1) = ep_set(n).xep(1) - 2*pi;
        ep_set_ext(m+6).xep(2) = ep_set(n).xep(2) - 4*pi;

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
%             xep = ep_set_ext(19).xep;
%             v = ep_set_ext(19).v;
%             perturb = 1e-2;
%             [~ , x_p] = ode45(@f_backward,[0,50],xep+v*perturb,odeset('RelTol',1e-5));
%             [~ , x_n] = ode45(@f_backward,[0,50],xep-v*perturb,odeset('RelTol',1e-5));
%             x_all = [flip(x_n,1);x_p];
%             plot(x_all(:,1),x_all(:,2),'b-','linewidth',1);
%             axis([0,2.5,0,3.5]);


for m = 1:length(ep_set)
    disp_v('Index',m);
    disp_v('Equilibrium',ep_set(m).xep);
    disp_v('Eigenvalue', ep_set(m).Lambda);
    disp_v('Eigenvector',ep_set(m).V);
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
function dfdt = f(x)  
        dfdt = f_reducedstate(x);
end
function dfdt = f_forward(t,x)
    dfdt = f(x);
end

function dfdt = f_backward(t,x)
    dfdt = -f(x);
end


