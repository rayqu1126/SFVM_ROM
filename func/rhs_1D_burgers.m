%% 1D Burgers' equations (no stochastic variable)
function dUdt = rhs_1D_burgers(~,U,params)
    N = params.N;
    dx = params.dx;
   
    % Reconstruct states in physical domain
    [UL,UR] = WENO_1D(U,N,dx);
    
    % Flux on the left interface
    Uim1  = [UR(N);UR(1:N-1)];
    Ui    = UL;
    Fim12 = (Uim1.^2 + Ui.^2)/4 - 0.5*max(abs(Ui),abs(Uim1)).*(Ui-Uim1);

    % Flux on the right interface (periodic)
    Fip12 = [Fim12(2:end); Fim12(1)];

    dUdt = -1/dx*(Fip12-Fim12);
end
