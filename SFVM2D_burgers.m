%% SFVM and ROM for 1D Burgers equations with 1D stochstic variables
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
ics = sin(2*pi*X)+ 0.5 * Y;

t0=0;
tmax=0.3;
tspan = [t0 tmax];

% Assemble parameters
dx = 1/Nx;
dy = 1/N;
params.N = N;
params.Nx = Nx;
params.dx = 1/Nx;
params.dy = 1/N;

% SFVM - state reconstruction
options = odeset('RelTol',1e-10,'AbsTol',1e-11);

tic
[~,U_state] = ode45(@(t,U) rhs_2D_burgers_state(t,U,params), tspan, ics, options);
toc

sol_state = reshape(U_state(end,:),Nx,N);

figure
surf(X,Y,sol_state)
xlabel('X'), ylabel('Y'), zlabel('SFVM(state)')


% SFVM - flux reconstruction

tic
[~,U_flux] = ode45(@(t,U) rhs_2D_burgers_flux(t,U,params), tspan, ics, options);
toc

sol_flux = reshape(U_flux(end,:),Nx,N);

figure
surf(X,Y,sol_flux)
xlabel('X'), ylabel('Y'), zlabel('SFVM(flux)')

er_flux = sum(abs(sol_state - sol_flux)...
    * (dx*dy),"all") / sum(abs(sol_state) * (dx*dy),"all");

disp("The relative L1 error (state and flux) is " + er_flux)

%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 200);
F_sample = zeros(2*N, Nx*length(tspan_sample));

quadL = y + 0.5*(-1/sqrt(3)) * dy;
quadR = y + 0.5*(1/sqrt(3)) * dy;   
quad_combined = reshape([quadL; quadR], [], 1)';
for i = 1:length(quad_combined)
    ics_sample = sin(2*pi*x) + 0.5 * quad_combined(i);
    F_sample(i,:)= sampling_burgers(Nx, ics_sample, tspan_sample);
end

% Construct reduced basis
Nmode = 30;
[Vn,s,~] = svd(F_sample);
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

[~,U_ROM] = ode45(@(t,U) rhs_2D_burgers_interp(t,U,params), tspan, ics, options);
%

sol_ROM = reshape(U_ROM(end,:),Nx,N);

figure
surf(X,Y,sol_ROM)
xlabel('X'), ylabel('Y'), zlabel('SFVM(ROM)')

er_ROM = sum(abs(sol_state - sol_ROM)...
    * (dx*dy),"all") / sum(abs(sol_state) * (dx*dy),"all");

disp("The relative L1 error (state and ROM) is " + er_ROM)