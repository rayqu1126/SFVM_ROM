%% Sampling of 1D Euler equation
function [F_snap] = sampling_euler(N, gamma, ics, tspan)
dx = 1/N;
params.dx = dx;
params.N = N;
params.gamma = gamma;

% Run 1D problem
options = odeset('RelTol',1e-6,'AbsTol',1e-8);
[~,U] = ode45(@(t,U) rhs_1D_euler(t,U,params), tspan, ics, options);

% Reconstruct states and collect flux
F_snap = [];
for i = 1:length(tspan)
    U_snap = U(i,:)';
    rho = U_snap(1:N);
    rhou = U_snap(N+1:2*N);
    E = U_snap(2*N+1,end);
    u = rhou ./ rho;
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);    

    [rhoL,rhoR] = WENO_1D_reflect(rho,N,dx,"rho");
    [uL,uR] = WENO_1D_reflect(u,N,dx,"u");
    [pL,pR] = WENO_1D_reflect(p,N,dx,"p");

    rhoim1 = [rhoR(1);rhoR(1:N-1)];
    rhoi = rhoL;
    uim1 = [-uR(1);uR(1:N-1)];
    ui = uL;
    pim1 = [pR(1);pR(1:N-1)];
    pi = pL;

    Eim1 = pim1 / (gamma - 1) + 0.5 * rhoim1 .* uim1.^2;
    Ei = pi / (gamma - 1) + 0.5 * rhoi .* ui.^2;    

    lambda = max(abs(ui)+sqrt(gamma*pi./rhoi),...
        abs(uim1)+sqrt(gamma*pim1./rhoim1));
    
    F_rho_im12 = (rhoi.*ui + rhoim1.*uim1)/2 - 0.5*lambda.*(rhoi-rhoim1);
    F_rhou_im12 = (rhoi.*ui.^2+pi + rhoim1.*uim1.^2+pim1)/2 - 0.5*lambda.*(rhoi .* ui- rhoim1 .* uim1);
    F_E_im12 = (ui.*(Ei + pi) + uim1.*(Eim1 + pim1))/2 - 0.5*lambda.*(Ei-Eim1);

    rhoip1 = [rhoL(2:N);rhoL(N)];
    rhoi = rhoR;
    uip1 = [uL(2:N);-uL(N)];
    ui = uR;
    pip1 = [pL(2:N);pL(N)];
    pi = pR;

    Eip1 = pip1 / (gamma - 1) + 0.5 * rhoip1 .* uip1.^2;
    Ei = pi / (gamma - 1) + 0.5 * rhoi .* ui.^2;    


    lambda = max(abs(ui)+sqrt(gamma*pi./rhoi),...
        abs(uip1)+sqrt(gamma*pip1./rhoip1));


    F_rho_ip12 = (rhoi.*ui + rhoip1.*uip1)/2 - 0.5*lambda.*(rhoip1-rhoi);
    F_rhou_ip12 = (rhoi.*ui.^2+pi + rhoip1.*uip1.^2+pip1)/2 - 0.5*lambda.*(rhoip1 .* uip1- rhoi .* ui);
    F_E_ip12 = (ui.*(Ei + pi) + uip1.*(Eip1 + pip1))/2 - 0.5*lambda.*(Eip1-Ei);

    Fim12 = [F_rho_im12; F_rhou_im12; F_E_im12];
    Fip12 = [F_rho_ip12; F_rhou_ip12; F_E_ip12];
    F = [Fim12;Fip12];
    F_snap = [F_snap; F];
end
F_snap = F_snap';
end
