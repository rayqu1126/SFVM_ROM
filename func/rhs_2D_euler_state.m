%% 1D Euler equations x 1D stochastic variable (state reconstruction)
function dUdt = rhs_2D_euler_state(~,U,params)
    Nx = params.Nx;
    N = params.N;
    dx = params.dx;
    gamma = params.gamma;

    % Unpack variables
    U = reshape(U,3*Nx,N);
    rho = U(1:Nx,:);
    rhou = U(Nx+1:2*Nx, :);
    E = U(2*Nx+1:end,:);
    
    u = rhou ./ rho;
    p = (gamma - 1) * (E - 0.5 * rho .* u.^2);
    
    % Reconstruct states in physical domain
    [rhoL,rhoR] = WENO_2D_reflect(rho,Nx,N,dx,"rho");
    [uL,uR] = WENO_2D_reflect(u,Nx,N,dx,"u");
    [pL,pR] = WENO_2D_reflect(p,Nx,N,dx,"p");

    % Reconstruct states in quadrature points
    a1 = 0.5*(-1/sqrt(3));
    a2 = 0.5*(1/sqrt(3));

    rhoL1 = (WENO_QP(N,a1) * rhoL')';
    rhoL2 = (WENO_QP(N,a2) * rhoL')';
    rhoR1 = (WENO_QP(N,a1) * rhoR')';
    rhoR2 = (WENO_QP(N,a2) * rhoR')';

    uL1 = (WENO_QP(N,a1) * uL')';
    uL2 = (WENO_QP(N,a2) * uL')';
    uR1 = (WENO_QP(N,a1) * uR')';
    uR2 = (WENO_QP(N,a2) * uR')';

    pL1 = (WENO_QP(N,a1) * pL')';
    pL2 = (WENO_QP(N,a2) * pL')';
    pR1 = (WENO_QP(N,a1) * pR')';
    pR2 = (WENO_QP(N,a2) * pR')';

   
    % Flux on the left interface (reflective wall)
    rhoim1_1  = [rhoR1(1,:);rhoR1(1:Nx-1,:)];
    rhoi_1    = rhoL1;
    rhoim1_2  = [rhoR2(1,:);rhoR2(1:Nx-1,:)];
    rhoi_2    = rhoL2;

    uim1_1  = [-uR1(1,:);uR1(1:Nx-1,:)];
    ui_1    = uL1;
    uim1_2  = [-uR2(1,:);uR2(1:Nx-1,:)];
    ui_2    = uL2;

    pim1_1  = [pR1(1,:);pR1(1:Nx-1,:)];
    pi_1    = pL1;
    pim1_2  = [pR2(1,:);pR2(1:Nx-1,:)];
    pi_2    = pL2;



    Eim1_1 = pim1_1 / (gamma - 1) + 0.5 * rhoim1_1 .* uim1_1.^2;
    Ei_1 = pi_1 / (gamma - 1) + 0.5 * rhoi_1 .* ui_1.^2; 


    Eim1_2 = pim1_2 / (gamma - 1) + 0.5 * rhoim1_2 .* uim1_2.^2;
    Ei_2 = pi_2 / (gamma - 1) + 0.5 * rhoi_2 .* ui_2.^2; 


    lambda1 = max(abs(ui_1)+real(sqrt(gamma*pi_1./rhoi_1)),...
        abs(uim1_1)+real(sqrt(gamma*pim1_1./rhoim1_1)));

    lambda2 = max(abs(ui_2)+real(sqrt(gamma*pi_2./rhoi_2)),...
        abs(uim1_2)+real(sqrt(gamma*pim1_2./rhoim1_2)));

    F_rho_im12_1 = (rhoi_1.*ui_1 + rhoim1_1.*uim1_1)/2 - 0.5*lambda1.*(rhoi_1-rhoim1_1);
    F_rho_im12_2 = (rhoi_2.*ui_2 + rhoim1_2.*uim1_2)/2 - 0.5*lambda2.*(rhoi_2-rhoim1_2);
    
    F_rhou_im12_1 = (rhoi_1.*ui_1.^2+pi_1 + rhoim1_1.*uim1_1.^2+pim1_1)/2 - 0.5*lambda1.*(rhoi_1.*ui_1-rhoim1_1.*uim1_1);
    F_rhou_im12_2 = (rhoi_2.*ui_2.^2+pi_2 + rhoim1_2.*uim1_2.^2+pim1_2)/2 - 0.5*lambda2.*(rhoi_2.*ui_2-rhoim1_2.*uim1_2);

    F_E_im12_1 = (ui_1.*(Ei_1 + pi_1) + uim1_1.*(Eim1_1 + pim1_1))/2 - 0.5*lambda1.*(Ei_1-Eim1_1);
    F_E_im12_2 = (ui_2.*(Ei_2 + pi_2) + uim1_2.*(Eim1_2 + pim1_2))/2 - 0.5*lambda2.*(Ei_2-Eim1_2);


    Fim12_1 = [F_rho_im12_1; F_rhou_im12_1; F_E_im12_1];
    Fim12_2 = [F_rho_im12_2; F_rhou_im12_2; F_E_im12_2];


    % Flux on the right interface (Reflective wall)
    rhoip1_1  = [rhoL1(2:Nx,:);rhoL1(Nx,:)];
    rhoi_1    = rhoR1;
    rhoip1_2  = [rhoL2(2:Nx,:);rhoL2(Nx,:)];
    rhoi_2    = rhoR2;

    uip1_1  = [uL1(2:Nx,:);-uL1(Nx,:)];
    ui_1    = uR1;
    uip1_2  = [uL2(2:Nx,:);-uL2(Nx,:)];
    ui_2    = uR2;

    pip1_1  = [pL1(2:Nx,:);pL1(Nx,:)];
    pi_1    = pR1;
    pip1_2  = [pL2(2:Nx,:);pL2(Nx,:)];
    pi_2    = pR2;

    
    Eip1_1 = pip1_1 / (gamma - 1) + 0.5 * rhoip1_1 .* uip1_1.^2;
    Ei_1 = pi_1 / (gamma - 1) + 0.5 * rhoi_1 .* ui_1.^2;    

    Eip1_2 = pip1_2 / (gamma - 1) + 0.5 * rhoip1_2 .* uip1_2.^2;
    Ei_2 = pi_2 / (gamma - 1) + 0.5 * rhoi_2 .* ui_2.^2;    

    lambda1 = max(abs(ui_1)+real(sqrt(gamma*pi_1./rhoi_1)),...
        abs(uip1_1)+real(sqrt(gamma*pip1_1./rhoip1_1)));

    lambda2 = max(abs(ui_2)+real(sqrt(gamma*pi_2./rhoi_2)),...
        abs(uip1_2)+real(sqrt(gamma*pip1_2./rhoip1_2)));
        
    Frhoip12_1 = (rhoi_1.*ui_1 + rhoip1_1.*uip1_1)/2 - 0.5*lambda1.*(rhoip1_1-rhoi_1);
    Frhoip12_2 = (rhoi_2.*ui_2 + rhoip1_2.*uip1_2)/2 - 0.5*lambda2.*(rhoip1_2-rhoi_2);
    
    Frhouip12_1 = (rhoi_1.*ui_1.^2+pi_1 + rhoip1_1.*uip1_1.^2+pip1_1)/2 - 0.5*lambda1.*(rhoip1_1.*uip1_1-rhoi_1.*ui_1);
    Frhouip12_2 = (rhoi_2.*ui_2.^2+pi_2 + rhoip1_2.*uip1_2.^2+pip1_2)/2 - 0.5*lambda2.*(rhoip1_2.*uip1_2-rhoi_2.*ui_2);

    FEip12_1 = (ui_1.*(Ei_1 + pi_1) + uip1_1.*(Eip1_1 + pip1_1))/2 - 0.5*lambda1.*(Eip1_1-Ei_1);
    FEip12_2 = (ui_2.*(Ei_2 + pi_2) + uip1_2.*(Eip1_2 + pip1_2))/2 - 0.5*lambda2.*(Eip1_2-Ei_2);


    Fip12_1 = [Frhoip12_1; Frhouip12_1; FEip12_1];
    Fip12_2 = [Frhoip12_2; Frhouip12_2; FEip12_2];




    dUdt = reshape(-1/dx*((Fip12_1+Fip12_2)/2-(Fim12_1+Fim12_2)/2), [],1);
end
