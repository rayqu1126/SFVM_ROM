%% 1D Burgers' equations x 1D stochastic variable (flux reconstruction)
function dUdt = rhs_2D_burgers_flux(~,U,params)
    Nx = params.Nx;
    N = params.N;
    dx = params.dx;

    % Reconstruct states in physical domain
    [UL,UR] = WENO_2D(reshape(U,Nx,N),Nx,N,dx);

    % Flux on the left interface
    Uim1  = [UR(Nx,:);UR(1:Nx-1,:)];
    Ui    = UL;
    Fim12 = (Uim1.^2 + Ui.^2)/4 - 0.5 * max(abs(Ui),abs(Uim1)).*(Ui-Uim1);
    % Flux on the right interface (periodic)
    Fip12 = [Fim12(2:end,:); Fim12(1,:)];

    % Reconstruct flux in quadrature points
    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));

    FyL1 = (WENO_QP(N,a1)*Fim12')';
    FyL2 = (WENO_QP(N,a2)*Fim12')';
    FyR1 = (WENO_QP(N,a1)*Fip12')';
    FyR2 = (WENO_QP(N,a2)*Fip12')';

    dUdt = reshape(-1/dx*((FyR1+FyR2)/2-(FyL1+FyL2)/2), [], 1);
end


