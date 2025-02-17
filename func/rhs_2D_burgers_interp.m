%% 1D Burgers' equations x 1D stochastic variable (interpolation method)
function dUdt = rhs_2D_burgers_interp(~,U,params)
    Nx = params.Nx;
    N = params.N;
    dx = params.dx;
    B = params.B;
    Vinv = params.Vinv;
    ids = params.ids;

    % Reconstruct states in physical domain
    [UL,UR] = WENO_2D(reshape(U,Nx,N),Nx,N,dx);
 
    % Flux on the left interface
    Uim1  = [UR(Nx,:);UR(1:Nx-1,:)];
    Ui    = UL;
    Fim12 = (Uim1.^2 + Ui.^2)/4 - 0.5*max(abs(Ui),abs(Uim1)).*(Ui-Uim1);

    % Reconstruct flux in quadrature points
    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));

    % Combine flux in order of quadrature points
    Fy1 = (WENO_QP(N,a1)*Fim12')';
    Fy2 = (WENO_QP(N,a2)*Fim12')';
    Fy_combined = zeros(2*N,Nx);
    for i = 1:Nx
        Fy_combined(:,i) =  reshape([Fy1(i,:); Fy2(i,:)], [], 1)';
    end

    % Flux integrals with reduced basis and HR indices
    intFim12 = B * Vinv * Fy_combined(ids,:); intFim12 = intFim12';
    intFip12 = [intFim12(2:end,:); intFim12(1,:)];

    dUdt = reshape( (-1/dx) * (intFip12 - intFim12), [], 1);
end




