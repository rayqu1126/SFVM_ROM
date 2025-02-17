%% WENO reconstruction at left and right edges for 1D array (periodic)
function [UL,UR] = WENO_1D(U,N,dx)
UL = zeros(N,1); UR = zeros(N,1);
dL = [1/3,2/3]; dR = [2/3,1/3];
    for i = 1:N
        if i == 1
            Ust1 = U(end,1);
            Ust2 = U(i,1);
            Ust3 = U(i+1,1);
        elseif i == N
            Ust1 = U(i-1,1);
            Ust2 = U(i,1);
            Ust3 = U(1,1);   
        else
            Ust1 = U(i-1,1);
            Ust2 = U(i,1);
            Ust3 = U(i+1,1);
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

        UL(i,1) = w0L * p0 + w1L * p1;

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

        UR(i,1) = w0L * p0 + w1L * p1;
    end

end