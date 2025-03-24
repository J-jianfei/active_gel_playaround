clear;
clc;
zeta = 33;    % Activity coefficient
alpha = 100;  % Nonlinear coefficient
gamma = 1;    % Gradient energy coefficient
kd = 200;      % Depolymerization rate
xi = 1;       % Friction coefficient
xi1 = 0.46;
% U = 10;      % Layer velocity
phi0 = 0.5;   % Base density

% pack parameters
params.zeta = zeta;
params.alpha = alpha;
params.gamma = gamma;
params.kd = kd;
params.xi = xi;
params.phi0 = phi0;
params.xi1 = xi1;



dt = 1e-5; dtsave = 0.005; j = 1;
tmax = 0.3;

dx = 10^-2;
L = 1;
Nx = round(L/dx);

x = linspace(0,L,Nx);

phi = linspace(0.5,0.5,Nx); %zeros(size(x)) + phi0;
phi(1) = phi(2);
phi(end) = phi(end-1);


dphidx = gradient(phi,dx);
d2phidx2 = gradient(dphidx,dx);
d3phidx3 = gradient(d2phidx2,dx);


U = -10; % -1/xi1 * (zeta*phi(end) - alpha*(phi(end)-phi0)^3+gamma*d2phidx2(end) - ...
    % zeta*phi(1) + alpha*(phi(1)-phi0)^3-gamma*d2phidx2(1));

params.U = U;


phi_prev(1,:) = phi;
phi_prev(2,:) = dphidx;
phi_prev(3,:) = d2phidx2;
phi_prev(4,:) = d3phidx3; 

params.dt = dt;

for t = 0:dt:tmax

    F{1} = griddedInterpolant(x,phi_prev(1,:));
    F{2} = griddedInterpolant(x,phi_prev(2,:));
    F{3} = griddedInterpolant(x,phi_prev(3,:));
    F{4} = griddedInterpolant(x,phi_prev(4,:));

    solinit = bvpinit(x,@(x)phi_guess(x,F));
    sol = bvp4c(@(x,phi) ode_system(x,phi,F,params),@(ya,yb) bcfun(ya,yb,params),solinit);

    phi = deval(sol,x,1);
    dphidx = deval(sol,x,2);
    d2phidx2 = deval(sol,x,3);
    d3phidx3 = deval(sol,x,4);

    if mod(t,dtsave) ==0
        phi_data(j,:) = phi;
        U_data(j) = U;
        j= j+1;
        disp(strcat("Current time"," ",num2str(t)));
        plot(x,phi);
        pause(0.05)
    end
   
    phi_prev(1,:) = phi;
    phi_prev(2,:) = dphidx;
    phi_prev(3,:) = d2phidx2;
    phi_prev(4,:) = d3phidx3;

    U = -1/xi1 * (zeta*phi(end) - alpha*(phi(end)-phi0)^3+gamma*d2phidx2(end) - ...
    zeta*phi(1) + alpha*(phi(1)-phi0)^3-gamma*d2phidx2(1));
    parms.U = U;



end


%%  functions

function phi_out = phi_guess(x,interplants)
    % phi_prev(1,:) = phi;
    % phi_prev(2,:) = dphi/dx;
    % phi_prev(3,:) = d2phi/dx2;
    % phi_prev(4,:) = d3phi/dx3;
    F  = interplants;
    phi_out = [
        F{1}(x);
        F{2}(x);
        F{3}(x);
        F{4}(x);
    ];
end

function dphidx = ode_system(x,phi,phi_prev_interp,params)

    zeta = params.zeta;
    alpha = params.alpha;
    gamma = params.gamma;
    kd = params.kd;
    xi = params.xi;
    U = params.U;
    phi0 = params.phi0;
    dt = params.dt;

    F = phi_prev_interp;

    phi_1 = phi(1); % phi
    phi_2 = phi(2); % phi'    
    phi_3 = phi(3); % phi''
    phi_4 = phi(4); % phi'''

    dphidx = zeros(4,1);
    dphidx(1) = phi_2;
    dphidx(2) = phi_3;
    dphidx(3) = phi_4;
    dphidx(4) = xi*U/gamma * phi_2 - zeta/gamma*phi_3 + 6*alpha/gamma*(phi_1-phi0)*phi_2^2 ...
        + 3*alpha/gamma*(phi_1-phi0)^2*phi_3 - xi*kd/gamma*(phi_1-phi0) - xi/gamma *(phi_1 - F{1}(x))/dt;

end

function res = bcfun(ya,yb,params)
    zeta = params.zeta;
    alpha = params.alpha;
    gamma = params.gamma;
    xi = params.xi;
    U = params.U;
    phi0 = params.phi0;

    res = [ya(2) ...
        ya(4) - (U*xi/gamma*ya(1)-(zeta/gamma - 3*alpha/gamma*(ya(1)-phi0)^2)*ya(2)) ...
        yb(2) ...
        yb(4) - (U*xi/gamma*yb(1)-(zeta/gamma - 3*alpha/gamma*(yb(1)-phi0)^2)*yb(2))];
end