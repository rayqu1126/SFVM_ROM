%% 1D Burgers' equations x 2D stochastic variable (flux reconstruction)
function dUdt = rhs_3D_burgers_flux(~,U,params)
    N = params.N;
    Nx = params.Nx;
    dx = params.dx;

    % Reconstruct states in physical domain
    U = reshape(U,Nx,N,N);
    UR = zeros(Nx,N,N); UL = zeros(Nx,N,N);
    for i = 1:N
        [UL(:,:,i),UR(:,:,i)] = WENO_2D(U(:,:,i),Nx,N,dx);
    end

    % Reconstruct flux on the left interface
    Uim1  = [UR(Nx,:,:);UR(1:Nx-1,:,:)];
    Ui    = UL;
    Fim12 = (Uim1.^2 + Ui.^2)/4 - 0.5*max(abs(Ui),abs(Uim1)).*(Ui-Uim1);

    % Reconstruct flux on the right interface
    Fip12 = [Fim12(2:end,:,:);Fim12(1,:,:)];

    % Reconstruct flux in quadrature points
    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));

    % Reconstruct flux in y 
    FyL1 = pagetranspose(pagemtimes(WENO_QP(N,a1),pagetranspose(Fim12)));
    FyL2 = pagetranspose(pagemtimes(WENO_QP(N,a2),pagetranspose(Fim12)));
    FyR1 = pagetranspose(pagemtimes(WENO_QP(N,a1),pagetranspose(Fip12)));
    FyR2 = pagetranspose(pagemtimes(WENO_QP(N,a2),pagetranspose(Fip12)));

    % Reconstruct flux in z
    FzL11 = permute(pagemtimes(WENO_QP(N,a1),permute(FyL1,[3 2 1])),[3 2 1]);
    FzL12 = permute(pagemtimes(WENO_QP(N,a2),permute(FyL1,[3 2 1])),[3 2 1]);
    FzL21 = permute(pagemtimes(WENO_QP(N,a1),permute(FyL2,[3 2 1])),[3 2 1]);
    FzL22 = permute(pagemtimes(WENO_QP(N,a2),permute(FyL2,[3 2 1])),[3 2 1]);
    FzR11 = permute(pagemtimes(WENO_QP(N,a1),permute(FyR1,[3 2 1])),[3 2 1]);
    FzR12 = permute(pagemtimes(WENO_QP(N,a2),permute(FyR1,[3 2 1])),[3 2 1]);
    FzR21 = permute(pagemtimes(WENO_QP(N,a1),permute(FyR2,[3 2 1])),[3 2 1]);
    FzR22 = permute(pagemtimes(WENO_QP(N,a2),permute(FyR2,[3 2 1])),[3 2 1]);


    dUdt = reshape(-1/dx*((FzR11+FzR12+FzR21+FzR22)/4 ...
        -(FzL11+FzL12+FzL21+FzL22)/4), [],1);
end


