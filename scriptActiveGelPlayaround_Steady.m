%%
%steady_state_solver
%%
%function steady_state_solver()
    % Parameters
    params.zeta = 30;    % Activity coefficient
    params.alpha = 100;  % Nonlinear coefficient
    params.gamma = 1;    % Gradient energy coefficient
    params.kd = 200;      % Depolymerization rate
    params.xi = 1;       % Friction coefficient
    params.U = 16.7;      % Layer velocity
    params.phi0 = 0.5;   % Base density

    % Domain
     L = 1;               % Length of the system
     xspan = linspace(0,L,10);     % Spatial domain
     solinit = bvpinit(xspan,[params.phi0;0;0;0]);

    % Initial guess for [phi, phi', phi'', phi''']
    % phi(0) = phi0; phi'(0) = 0; no-flux-derived estimate for phi''(0); phi'''(0) = 0
    % y0 = [params.phi0; 0; (params.zeta / params.gamma) * params.phi0; 0];

    % Solve using ode45
    sol = bvp5c(@(x,y) ode_system(x,y,params),@(ya,yb) bcfun(ya,yb,params),solinit);
    % Plot results
    
    figure;
    plot(sol.x, sol.y(1,:), 'b-', 'LineWidth', 2);
    xlabel('x');
    ylabel('\phi');
    title('Steady-State Solution for \phi(x)');
    grid on;
%end

function dydx = ode_system(x, y, params)
    % Extract parameters
    zeta = params.zeta;
    alpha = params.alpha;
    gamma = params.gamma;
    kd = params.kd;
    xi = params.xi;
    U = params.U;
    phi0 = params.phi0;

    % Variables
    y1 = y(1); % phi
    y2 = y(2); % phi'
    y3 = y(3); % phi''
    y4 = y(4); % phi'''

    % System of ODEs
    dydx = zeros(4, 1);
    dydx(1) = y2; % phi' = y2
    dydx(2) = y3; % phi'' = y3
    dydx(3) = y4; % phi''' = y4
    dydx(4) = xi*U/gamma * y2 - zeta/gamma*y3+6*alpha/gamma*(y1-phi0)*y2^2 ...
        + 3*alpha/gamma*(y1-phi0)^2*y3 - xi*kd/gamma*(y1-phi0);
end

function res = bcfun(ya,yb,params)
    zeta = params.zeta;
    alpha = params.alpha;
    gamma = params.gamma;
    U = params.U;
    phi0 = params.phi0;

    res = [ya(2) ...
        ya(4) - (U/gamma*ya(1)-(zeta/gamma - 3*alpha/gamma*(ya(1)-phi0)^2)*ya(2)) ...
        yb(2) ...
        yb(4) - (U/gamma*yb(1)-(zeta/gamma - 3*alpha/gamma*(yb(1)-phi0)^2)*yb(2))];

end