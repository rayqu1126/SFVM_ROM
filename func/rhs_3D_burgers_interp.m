%% 1D Burgers' equations x 2D stochastic variable (interpolation method)
function dUdt = rhs_3D_burgers_interp(~,U,params)
  N = params.N;
  Nx = params.Nx;
  dx = params.dx;
  B = params.B;
  Vinv = params.Vinv;
  ids = params.ids;

  % Reconstruct states in physical domain
  U = reshape(U,Nx,N,N);
  UR = zeros(Nx,N,N); UL = zeros(Nx,N,N);
  for i = 1:N
      [UL(:,:,i),UR(:,:,i)] = WENO_2D(U(:,:,i),Nx,N,dx);
  end
  
  Uim1  = [UR(Nx,:,:);UR(1:Nx-1,:,:)];
  Ui    = UL;
  Fim12 = (Uim1.^2 + Ui.^2)/4 - 0.5*max(abs(Ui),abs(Uim1)).*(Ui-Uim1);

  % Flux in qudrature points
  a1 = 0.5*(-1/sqrt(3));
  a2 = 0.5*(1/sqrt(3));

  FyL1 = pagetranspose(pagemtimes(WENO_QP(N,a1),pagetranspose(Fim12)));
  FyL2 = pagetranspose(pagemtimes(WENO_QP(N,a2),pagetranspose(Fim12)));

  FzL11 = permute(pagemtimes(WENO_QP(N,a1),permute(FyL1,[3 2 1])),[3 2 1]);
  FzL12 = permute(pagemtimes(WENO_QP(N,a2),permute(FyL1,[3 2 1])),[3 2 1]);
  FzL21 = permute(pagemtimes(WENO_QP(N,a1),permute(FyL2,[3 2 1])),[3 2 1]);
  FzL22 = permute(pagemtimes(WENO_QP(N,a2),permute(FyL2,[3 2 1])),[3 2 1]);

  % Combine flux at quadrature points (with same arrangement as sampling)

  FzCombined = zeros(4*N^2,Nx);
  for i = 1:Nx
    Fz11 = reshape(FzL11(i,:,:),N,N);
    Fz12 = reshape(FzL12(i,:,:),N,N);
    Fz21 = reshape(FzL21(i,:,:),N,N);
    Fz22 = reshape(FzL22(i,:,:),N,N);
    FzCombined(:,i) =  reshape([reshape(Fz11,1,N^2); reshape(Fz21,1,N^2); ...
                    reshape(Fz12,1,N^2); reshape(Fz22,1,N^2)],4*N^2,1);
  end
  
  % Flux integrals with reduced basis and HR indices
  intFim12 = B * Vinv * FzCombined(ids,:); 
  intFim12 = intFim12';
  intFim12 = reshape(intFim12,Nx,N,N);


  intFip12 = [intFim12(2:end,:,:);intFim12(1,:,:)];

  dUdt = reshape(-1/dx*(intFip12 - intFim12), [],1);
end

