%% SFVM and ROM for 1D Burgers equations with 2D stochstic variables
clear
close all
addpath('func')

% Set default plot
set(0, 'defaultaxesfontsize',24,'defaultaxeslinewidth',4,...
       'defaultlinelinewidth',4,'defaultpatchlinewidth',4,...
       'defaulttextfontsize',24,'defaulttextinterpreter','latex');
set(0, 'DefaultFigurePosition',[100 100 1100 600]);

% Initialize mesh (Nx: physical N: stochastic)
Nx = 64; 
N = 8;

xpts = linspace(0,1,Nx+1);
ypts = linspace(0,1,N+1);
zpts = linspace(0,1,N+1);

x = 0.5*(xpts(1:Nx)+xpts(2:Nx+1));
y = 0.5*(ypts(1:N)+ypts(2:N+1));
z = 0.5*(zpts(1:N)+zpts(2:N+1));

[Y,X,Z] = meshgrid(y,x,z);

% Initialize solution array
ics = 0.1 * Z + sin(2*pi*X) + 0.2 * sin(2*pi*Y);

t0=0;
tmax=0.2;
tspan = [t0 tmax];

dx = 1/Nx;
dy = 1/N;
dz = 1/N;

% Assemble parameters
params.N = N;
params.Nx = Nx;
params.dx = dx;
params.dy = dy;
params.dz = dz;

% SFVM - state reconstruction
options = odeset('RelTol',1e-10,'AbsTol',1e-11);

tic
[~,U_state] = ode45(@(t,U) rhs_3D_burgers_state(t,U,params), tspan, ics, options);
toc

sol_state =reshape(U_state(end,:),Nx,N,N);

[YY,XX] = meshgrid(y,x);
figure
surf(XX,YY,sol_state(:,:,end))
xlabel('X'), ylabel('Y'), zlabel('SFVM (state)')

% SFVM - flux reconstruction
[~,U_flux] = ode45(@(t,U) rhs_3D_burgers_flux(t,U,params), tspan, ics, options);

sol_flux=reshape(U_flux(end,:),Nx,N,N);

figure
surf(XX,YY,sol_flux(:,:,end))
xlabel('X'), ylabel('Y'), zlabel('SFVM (flux)')

er_flux = sum(abs(sol_state - sol_flux)...
    * (dx*dy*dz),"all") / sum(abs(sol_state) * (dx*dy*dz),"all");

disp("The relative L1 error (state and flux) is " + er_flux)

%% ROM - POD
% 1D sampling
tspan_sample = linspace(t0, tmax, 200);
F_sample = zeros(4*N^2,Nx*length(tspan_sample));

% Coordinates y and z of quad points volume by volume
yquadL = y + 0.5*(-1/sqrt(3)) * dy;
yquadR = y + 0.5*(1/sqrt(3)) * dy;   
yquad = reshape([yquadL; yquadR; yquadL; yquadR], [], 1)';
yquad = repmat(yquad,1,N);

zquadL = z + 0.5*(-1/sqrt(3)) * dz;
zquadR = z + 0.5*(1/sqrt(3)) * dz;
zquad = repmat([zquadL; zquadL; zquadR; zquadR],N,1);
zquad = reshape(zquad,[],1)';

for i = 1:length(yquad)
    ics_sample = 0.1 * zquad(i) + sin(2*pi*x') + ...
        0.2 * sin(2*pi*yquad(i));
    F_sample(i,:) = sampling_burgers(Nx, ics_sample, tspan_sample);
end

% Construct reduced basis
Nmode = 30;
[Vn,s,~] = svd(F_sample);
Vn = Vn(:,1:Nmode);

% Integrals of bases
Bn = zeros(N^2,Nmode);
for i = 1:N^2
    Bn(i,:) = (Vn(4*i-3,:) +Vn(4*i-2,:)+Vn(4*i-1,:)+ Vn(4*i,:)) /4;
end

% Q-DEIM for hyper-reduction
[~,~,P] = qr(Vn',"vector");
NHR = Nmode;
ids = P(1:NHR);
Vn_HR = Vn(ids,:);

%
params.Vinv = pinv(Vn_HR);
params.B = Bn;
params.ids = ids;

tic
[~,U_ROM] = ode45(@(t,U) rhs_3D_burgers_interp(t,U,params), tspan, ics, options);
toc

sol_ROM = reshape(U_ROM(end,:),Nx,N,N);

er_ROM = sum(abs(sol_ROM - sol_state)...
    * (dx*dy*dz),"all") / sum(abs(sol_state) * (dx*dy*dz),"all");

disp("The relative L1 error (state and ROM) is " + er_ROM)

figure
surf(XX,YY,sol_flux(:,:,end))
xlabel('X'), ylabel('Y'), zlabel('SFVM (flux)')
% 

% Compute solution mean and std
sol_mean = sum(sum(sol_ROM, 2), 3) * dy * dz;  
sol_diff = sol_ROM - sol_mean;  
sol_var = sum(sum(sol_diff.^2, 2), 3) / (N^2 - 1); 
sol_std = sqrt(sol_var); 
% 
figure
plot(x, sol_mean, LineWidth= 2, color = "r")
hold on 
plot(x,sol_mean + sol_std, LineWidth= 2, color = "r", LineStyle= "--")
hold on
plot(x,sol_mean - sol_std, LineWidth=2, color = "r", LineStyle= "--")
title(['Nmode = ' num2str(Nmode) ', HR nodes = ' num2str(NHR)])
title('Mean and std of ROM solution')


