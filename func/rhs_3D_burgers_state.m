%% 1D Burgers' equations x 2D stochastic variable (state reconstruction)
function dUdt = rhs_3D_burgers_state(~,U,params)
    N = params.N;
    Nx = params.Nx;
    dx = params.dx;

    % Reconstruct states in physical domain
    U = reshape(U,Nx,N,N);
    UR = zeros(Nx,N,N); UL = zeros(Nx,N,N);
    for i = 1:N
        [UL(:,:,i),UR(:,:,i)] = WENO_2D(U(:,:,i),Nx,N,dx);
    end

    % Reconstruct states in quadrature points
    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));

    % Reconstruct states in y
    UyL1 = pagetranspose(pagemtimes(WENO_QP(N,a1),pagetranspose(UL)));
    UyL2 = pagetranspose(pagemtimes(WENO_QP(N,a2),pagetranspose(UL)));
    UyR1 = pagetranspose(pagemtimes(WENO_QP(N,a1),pagetranspose(UR)));
    UyR2 = pagetranspose(pagemtimes(WENO_QP(N,a2),pagetranspose(UR)));

    % Reconstruct states in z
    UzL11 = permute(pagemtimes(WENO_QP(N,a1),permute(UyL1,[3,2,1])),[3,2,1]);
    UzL12 = permute(pagemtimes(WENO_QP(N,a2),permute(UyL1,[3,2,1])),[3,2,1]);
    UzL21 = permute(pagemtimes(WENO_QP(N,a1),permute(UyL2,[3,2,1])),[3,2,1]);
    UzL22 = permute(pagemtimes(WENO_QP(N,a2),permute(UyL2,[3,2,1])),[3,2,1]);
    UzR11 = permute(pagemtimes(WENO_QP(N,a1),permute(UyR1,[3,2,1])),[3,2,1]);
    UzR12 = permute(pagemtimes(WENO_QP(N,a2),permute(UyR1,[3,2,1])),[3,2,1]);
    UzR21 = permute(pagemtimes(WENO_QP(N,a1),permute(UyR2,[3,2,1])),[3,2,1]);
    UzR22 = permute(pagemtimes(WENO_QP(N,a2),permute(UyR2,[3,2,1])),[3,2,1]);



    % Flux on the left interface
    Uim1_1  = [UzR11(Nx,:,:);UzR11(1:Nx-1,:,:)];
    Ui_1    = UzL11;
    Fim12_1 = (Uim1_1.^2 + Ui_1.^2)/4 - 0.5*max(abs(Ui_1),abs(Uim1_1))...
        .*(Ui_1-Uim1_1);


    Uim1_2  = [UzR12(Nx,:,:);UzR12(1:Nx-1,:,:)];
    Ui_2    = UzL12;
    Fim12_2 = (Uim1_2.^2 + Ui_2.^2)/4 - 0.5*max(abs(Ui_2),abs(Uim1_2))...
        .*(Ui_2-Uim1_2);

    Uim1_3  = [UzR21(Nx,:,:);UzR21(1:Nx-1,:,:)];
    Ui_3    = UzL21;
    Fim12_3 = (Uim1_3.^2 + Ui_3.^2)/4 - 0.5*max(abs(Ui_3),abs(Uim1_3))...
        .*(Ui_3-Uim1_3);

    Uim1_4  = [UzR22(Nx,:,:);UzR22(1:Nx-1,:,:)];
    Ui_4    = UzL22;
    Fim12_4 = (Uim1_4.^2 + Ui_4.^2)/4 - 0.5*max(abs(Ui_4),abs(Uim1_4))...
        .*(Ui_4-Uim1_4);

    % Flux on the right interface (periodic)
    Fip12_1 = [Fim12_1(2:end,:,:);Fim12_1(1,:,:)];
    Fip12_2 = [Fim12_2(2:end,:,:);Fim12_2(1,:,:)];
    Fip12_3 = [Fim12_3(2:end,:,:);Fim12_3(1,:,:)];
    Fip12_4 = [Fim12_4(2:end,:,:);Fim12_4(1,:,:)];

    
    dUdt = reshape(-1/dx*((Fip12_1+Fip12_2+Fip12_3+Fip12_4)/4-...
        (Fim12_1+Fim12_2+Fim12_3+Fim12_4)/4), [],1);
end


