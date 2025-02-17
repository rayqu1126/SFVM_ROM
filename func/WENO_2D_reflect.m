%% WENO reconstruction at left and right edges for 2D array (reflective wall)
function [UL,UR] = WENO_2D_reflect(U,Nx,N,dx,var_type)
% (assume x is physical domain)
UL = zeros(Nx,N); UR = zeros(Nx,N);
dL = [1/3,2/3]; dR = [2/3,1/3];
for j = 1:N
    for i = 1:Nx
        if i == 1
            Ust1 = (var_type == "u") * (-U(1,j)) + (var_type ~= "u") * U(1,j); 
            Ust2 = U(i,j);
            Ust3 = U(i+1,j);
        elseif i == Nx
            Ust1 = U(i-1,j);
            Ust2 = U(i,j);
            Ust3 = (var_type == "u") * (-U(Nx,j)) + (var_type ~= "u") * U(Nx,j);   
        else
            Ust1 = U(i-1,j);
            Ust2 = U(i,j);
            Ust3 = U(i+1,j);
        end
        
        der0 = (Ust3 - Ust2) / dx;
        der1 = (Ust2 - Ust1) / dx;


        % smoothness indicator
        is0 = (Ust3 - Ust2)^2;
        is1 = (Ust2 - Ust1)^2;
        
        % reconstruct at left
        a0 = dL(1) / (eps + is0)^2;
        a1 = dL(2) / (eps + is1)^2;
        w0L = a0 / (a0 + a1);
        w1L = a1 / (a0 + a1);

        p0 = Ust2 - 0.5*der0*dx;
        p1 = Ust2 - 0.5*der1*dx;

        UL(i,j) = w0L * p0 + w1L * p1;

        % reconstruct at right
        is0 = (Ust3 - Ust2)^2;
        is1 = (Ust2 - Ust1)^2;

        % reconstruct at left
        a0 = dR(1) / (eps + is0)^2;
        a1 = dR(2) / (eps + is1)^2;
        w0L = a0 / (a0 + a1);
        w1L = a1 / (a0 + a1);

        p0 = Ust2 + 0.5*der0*dx;
        p1 = Ust2 + 0.5*der1*dx;

        UR(i,j) = w0L * p0 + w1L * p1;
    end
end

end