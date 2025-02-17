%% Sampling of 1D Burgers' equation
function [F_snap] = sampling_burgers(N, ics, tspan)
params.dx = 1/N;
params.N = N;

% Run 1D problem
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
[~,U] = ode45(@(t,U) rhs_1D_burgers(t,U,params), tspan, ics, options);

% Reconstruct states and collect flux
Fim12 = [];
for i = 1:length(tspan)
    [UL,UR] = WENO_1D(U(i,:)', N, params.dx);
    Uim1  = [UR(N,:);UR(1:N-1,:)]; 
    Ui    = UL;
    Fim12 = [Fim12; (Uim1.^2 + Ui.^2)/4 - ...
        0.5 * max(abs(Ui),abs(Uim1)).*(Ui-Uim1)];
end
F_snap = Fim12';
end
