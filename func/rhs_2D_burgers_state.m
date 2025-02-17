%% 1D Burgers' equations x 1D stochastic variable (state reconstruction)
function dUdt = rhs_2D_burgers_state(~,U,params)
    Nx = params.Nx;
    N = params.N;
    dx = params.dx;

    % Reconstruct states in physical domain
    [UL,UR] = WENO_2D(reshape(U,Nx,N),Nx,N,dx);

    % Reconstruct states in quadrature points
    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));

    UyL1 = (WENO_QP(N,a1) * UL')';
    UyL2 = (WENO_QP(N,a2) * UL')';
    UyR1 = (WENO_QP(N,a1) * UR')';
    UyR2 = (WENO_QP(N,a2) * UR')';

   
    % Flux on the left and right interface (periodic)
    Uim1_1  = [UyR1(Nx,:);UyR1(1:Nx-1,:)];
    Ui_1    = UyL1;
    Fim12_1 = (Uim1_1.^2 + Ui_1.^2)/4 - 0.5*max(abs(Ui_1),abs(Uim1_1)).*(Ui_1-Uim1_1);

    Uim1_2  = [UyR2(Nx,:);UyR2(1:Nx-1,:)];
    Ui_2    = UyL2;
    Fim12_2 = (Uim1_2.^2 + Ui_2.^2)/4 - 0.5*max(abs(Ui_2),abs(Uim1_2)).*(Ui_2-Uim1_2);


    %     
    Fip12_1 = [Fim12_1(2:end,:);Fim12_1(1,:)];
    Fip12_2 = [Fim12_2(2:end,:);Fim12_2(1,:)];


    dUdt = reshape(-1/dx*((Fip12_1+Fip12_2)/2-(Fim12_1+Fim12_2)/2), [], 1);
end