%% 1D Euler equations x 1D stochastic variable (interpolation method)
function dUdt = rhs_2D_euler_interp(~,U,params)
  N = params.N;
  Nx = params.Nx;
  dx = params.dx;
  B = params.B;
  Vinv = params.Vinv;
  ids = params.ids;
  gamma = params.gamma;


  % Unpack variables
  U = reshape(U,3 * Nx,N);
  rho = U(1:Nx,:);
  rhou = U(Nx+1:2*Nx, :);
  E = U(2*Nx+1:end,:);
    
  u = rhou ./ rho;
  p = (gamma - 1) * (E - 0.5 * rho .* u.^2);

  % Reconstruct states in physical domain
  [rhoL,rhoR] = WENO_2D_reflect(rho,Nx,N,dx,"rho");
  [uL,uR] = WENO_2D_reflect(u,Nx,N,dx,"u");
  [pL,pR] = WENO_2D_reflect(p,Nx,N,dx,"p");

  % Flux on left interface (reflective wall)
  rhoim1 = [rhoR(1,:);rhoR(1:Nx-1,:)];
  rhoi = rhoL;
  uim1 = [-uR(1,:);uR(1:Nx-1,:)];
  ui = uL;
  pim1 = [pR(1,:);pR(1:Nx-1,:)];
  pi = pL;

  Eim1 = pim1 / (gamma - 1) + 0.5 * rhoim1 .* uim1.^2;
  Ei = pi / (gamma - 1) + 0.5 * rhoi .* ui.^2;    

  lambda = max(abs(ui)+sqrt(gamma*pi./rhoi),...
        abs(uim1)+sqrt(gamma*pim1./rhoim1));
    
  F_rho_im12 = (rhoi.*ui + rhoim1.*uim1)/2 - 0.5*lambda.*(rhoi-rhoim1);
  F_rhou_im12 = (rhoi.*ui.^2+pi + rhoim1.*uim1.^2+pim1)/2 - 0.5*lambda.*(rhoi .* ui- rhoim1 .* uim1);
  F_E_im12 = (ui.*(Ei + pi) + uim1.*(Eim1 + pim1))/2 - 0.5*lambda.*(Ei-Eim1);

  % Flux on right interface (reflective wall)
  rhoip1 = [rhoL(2:Nx,:);rhoL(Nx,:)];
  rhoi = rhoR;
  uip1 = [uL(2:Nx,:);-uL(Nx,:)];
  ui = uR;
  pip1 = [pL(2:Nx,:);pL(Nx,:)];
  pi = pR;

  Eip1 = pip1 / (gamma - 1) + 0.5 * rhoip1 .* uip1.^2;
  Ei = pi / (gamma - 1) + 0.5 * rhoi .* ui.^2;    

   lambda = max(abs(ui)+sqrt(gamma*pi./rhoi),...
        abs(uip1)+sqrt(gamma*pip1./rhoip1));


  F_rho_ip12 = (rhoi.*ui + rhoip1.*uip1)/2 - 0.5*lambda.*(rhoip1-rhoi);
  F_rhou_ip12 = (rhoi.*ui.^2+pi + rhoip1.*uip1.^2+pip1)/2 - 0.5*lambda.*(rhoip1 .* uip1- rhoi .* ui);
  F_E_ip12 = (ui.*(Ei + pi) + uip1.*(Eip1 + pip1))/2 - 0.5*lambda.*(Eip1-Ei);
   
  % Reconstruct flux in quadrature points
  a1 = 0.5*(-1/sqrt(3)); 
  a2 = 0.5*(1/sqrt(3));

  F_rho_im12_1 = (WENO_QP(N,a1)*F_rho_im12')';
  F_rho_im12_2 = (WENO_QP(N,a2)*F_rho_im12')';
  F_rhou_im12_1 = (WENO_QP(N,a1)*F_rhou_im12')';
  F_rhou_im12_2 = (WENO_QP(N,a2)*F_rhou_im12')';
  F_E_im12_1 = (WENO_QP(N,a1)*F_E_im12')';
  F_E_im12_2 = (WENO_QP(N,a2)*F_E_im12')';

  Fim12_1 = [F_rho_im12_1; F_rhou_im12_1; F_E_im12_1];
  Fim12_2 = [F_rho_im12_2; F_rhou_im12_2; F_E_im12_2];

  Fim12_combined = zeros(2*N,3*Nx);
  for i = 1:3*Nx
      Fim12_combined(:,i) =  reshape([Fim12_1(i,:); Fim12_2(i,:)], [], 1)';
  end


  F_rho_ip12_1 = (WENO_QP(N,a1)*F_rho_ip12')';
  F_rho_ip12_2 = (WENO_QP(N,a2)*F_rho_ip12')';
  F_rhou_ip12_1 = (WENO_QP(N,a1)*F_rhou_ip12')';
  F_rhou_ip12_2 = (WENO_QP(N,a2)*F_rhou_ip12')';
  F_E_ip12_1 = (WENO_QP(N,a1)*F_E_ip12')';
  F_E_ip12_2 = (WENO_QP(N,a2)*F_E_ip12')';

  Fip12_1 = [F_rho_ip12_1; F_rhou_ip12_1; F_E_ip12_1];
  Fip12_2 = [F_rho_ip12_2; F_rhou_ip12_2; F_E_ip12_2];
    
  % Combine flux in order of quadrature
  Fip12_combined = zeros(2*N,3*Nx);
  for i = 1:3*Nx
    Fip12_combined(:,i) =  reshape([Fip12_1(i,:); Fip12_2(i,:)], [], 1)';
  end

  % Flux integrals with reduced basis and HR indices
  intFim12 = B * Vinv * Fim12_combined(ids,:); intFim12 = intFim12';
  intFip12 = B * Vinv * Fip12_combined(ids,:); intFip12 = intFip12';

  dUdt = reshape(-1/dx*(intFip12 - intFim12), [],1);
end