%% WENO reconstruction at quadrature ppints (outflow BCs)
function [Rec_a] = WENO_QP(N,a)
% Ghost cells are absorbed into the first and last cells.
d0 = 0.5*(a*a + a - 1./12.)/a;
d1 = 1.0 - d0;

Rec_a = zeros(N,N);

for i = 1:N
    if (i > 1)
        Rec_a(i,i-1) = -a*d1; 
    else
        Rec_a(i,i) = Rec_a(i,i) - a*d1; 
    end

    Rec_a(i,i) = Rec_a(i,i) + d0*(1.-a) + d1*(1.+a);

    if (i < N)
        Rec_a(i,i+1) = a*d0; 
    else
        Rec_a(i,i) = Rec_a(i,i) + a*d0; 
    end
end
end

