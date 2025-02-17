%% WENO reconstruction at left and right edges for 1D array (reflective wall)
function [UL,UR] = WENO_1D_reflect(U,N,dx,var_type)
UL = zeros(N,1); UR = zeros(N,1);
dL = [1/3,2/3]; dR = [2/3,1/3];
    for i = 1:N
        % Implement Reflective Wall Conditions for Ghost Cells
        if i == 1
            Ust1 = (var_type == "u") * (-U(1)) + (var_type ~= "u") * U(1); % Reflect velocity
            Ust2 = U(i);
            Ust3 = U(i+1);
        elseif i == N
            Ust1 = U(i-1);
            Ust2 = U(i);
            Ust3 = (var_type == "u") * (-U(N)) + (var_type ~= "u") * U(N); % Reflect velocity
        else
            Ust1 = U(i-1);
            Ust2 = U(i);
            Ust3 = U(i+1);
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

        UL(i) = w0L * p0 + w1L * p1;

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

        UR(i) = w0L * p0 + w1L * p1;
    end

end
