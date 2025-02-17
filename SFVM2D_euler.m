%% SFVM and ROM for 1D Euler equations with 1D stochastic variables
clear
close all
addpath('func')

% Set default plot
set(0, 'defaultaxesfontsize',24,'defaultaxeslinewidth',4,...
       'defaultlinelinewidth',4,'defaultpatchlinewidth',4,...
       'defaulttextfontsize',24,'defaulttextinterpreter','latex');
set(0, 'DefaultFigurePosition',[100 100 1100 600]);

% Initialize mesh (Nx: physical N: stochastic)
Nx = 128;
N = 16;

xpts = linspace(0,1,Nx+1);
ypts = linspace(0,1,N+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:N)+ypts(2:N+1));

[Y,X] = meshgrid(y,x);

% Initialize solution array
t0=0;
tmax=0.2;
tspan = [t0 tmax];

u = zeros(Nx,N);
% rho = ones(Nx,N) * 0.125;
% p = ones(Nx,N) * 0.1;
rho = ones(Nx,N) * 0.15;
p = ones(Nx,N) * 0.15;

rho(X < 0.5 + 0.1*Y) = 1; 
p(X < 0.5 + 0.1*Y) = 1;   
gamma = 1.4;
E = p / (gamma - 1) + 0.5 * rho .* u.^2;

ics = [rho; rho .* u; E];

% Assemble parameters
dx = 1/Nx;
dy = 1/N;

params.N = N;
params.Nx = Nx;
params.dx = dx;
params.dy = dy;
params.gamma = gamma;

% SFVM - state reconstruction 
options = odeset('NonNegative', numel(ics), 'RelTol',1e-8,'AbsTol',1e-10);
tic
[~,U_state] = ode45(@(t,U) rhs_2D_euler_state(t,U,params), tspan, ics, options);
toc

sol_state = reshape(U_state(end,:),3*Nx,N);
sol_state_rho = sol_state(1:Nx,:);


figure
surf(X,Y,sol_state_rho)
xlabel('X'), ylabel('Y'), zlabel('SFVM(state)')

% SFVM - flux reconstruction
tic
[~,U_flux] = ode45(@(t,U) rhs_2D_euler_flux(t,U,params), tspan, ics, options);
toc

sol_flux = reshape(U_flux(end,:),3*Nx,N);
sol_flux_rho = sol_flux(1:Nx,:);


figure
surf(X,Y,sol_flux_rho)
xlabel('X'), ylabel('Y'), zlabel('SFVM(flux)')


er_flux = sum(abs(sol_state - sol_flux) * (dx*dy),"all") ...
    / sum(abs(sol_state) * (dx*dy),"all");

disp("The relative L1 error (state and flux) is " + er_flux)

%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 100);
quadL = y + 0.5*(-1/sqrt(3)) * dy;
quadR = y + 0.5*(1/sqrt(3)) * dy;  

F_sample = zeros(2*N, 2*3*Nx*length(tspan_sample));
quad_combined = reshape([quadL; quadR], [], 1)';
for i = 1:length(quad_combined)
    u = zeros(Nx,1);
    rho = ones(Nx,1) * 0.15;
    p = ones(Nx,1) * 0.15;
    rho(x < 0.5+0.5*quad_combined(i)) = 1; 
    p(x < 0.5+0.5*quad_combined(i)) = 1;   
    E = p / (gamma - 1) + 0.5 * rho .* u.^2;    
    ics_sample = [rho; rho.*u; E];
    F_sample(i,:)= sampling_euler(Nx, gamma, ics_sample, tspan_sample);
end

% Construct reduced basis
Nmode = 30;
[Vn,s,~] = svd(F_sample,"econ"); % ECON-SVD for large matrix
Vn = Vn(:,1:Nmode);

% Integral of basis
Bn = zeros(N,Nmode);
for i = 1:N
    Bn(i,:) = (Vn(2*i-1,:) + Vn(2*i,:)) /2;
end

% Q-DEIM hyper-reduction
[~,~,P] = qr(Vn',"vector");

NHR = Nmode;
ids = P(1:NHR);

Vn_HR = Vn(ids,:);

%
params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;

[~,U_ROM] = ode45(@(t,U) rhs_2D_euler_interp(t,U,params), tspan, ics, options);
%

sol_ROM = reshape(U_ROM(end,:),3*Nx,N);
sol_ROM_rho = sol_ROM(1:Nx,:);

figure
surf(X,Y,sol_ROM_rho)
xlabel('X'), ylabel('Y'), zlabel('SFVM(ROM)')

er_ROM = sum(abs(sol_state - sol_ROM) * (dx*dy),"all") ...
    / sum(abs(sol_state) * (dx*dy),"all");

disp("The relative L1 error (state and ROM) is " + er_ROM)