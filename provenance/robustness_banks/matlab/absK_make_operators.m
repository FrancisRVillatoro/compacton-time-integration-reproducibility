function ops = absK_make_operators(N,dx,Lmax)
% Sparse periodic Padé operators.
% AA: common denominator.
% BB/AA: first derivative.
% CC/AA: third derivative.
% DD/AA: fourth derivative.
% SS/AA: second derivative, retained only for legacy diagnostics.

unos = ones(N,1);
ceros = zeros(N,1);

AA = 1/120*(spdiags([unos, 26*unos, 66*unos, 26*unos, unos],-2:2,N,N) + ...
            spdiags([unos, 26*unos, unos, 26*unos],[N-2,N-1,-N+2,-N+1],N,N));

BB = 1/(24*dx)*(spdiags([-unos, -10*unos, ceros, 10*unos, unos],-2:2,N,N) + ...
                spdiags([-unos, -10*unos, unos, 10*unos],[N-2,N-1,-N+2,-N+1],N,N));

CC = 1/(2*dx^3)*(spdiags([-unos, 2*unos, ceros, -2*unos, unos],-2:2,N,N) + ...
                 spdiags([-unos, 2*unos, unos, -2*unos],[N-2,N-1,-N+2,-N+1],N,N));

SS = 1/(6*dx^2)*(spdiags([unos, 2*unos, -6*unos, 2*unos, unos],-2:2,N,N) + ...
                 spdiags([unos, 2*unos, unos, 2*unos],[N-2,N-1,-N+2,-N+1],N,N));

DD = 1/dx^4*(spdiags([unos, -4*unos, 6*unos, -4*unos, unos],-2:2,N,N) + ...
             spdiags([unos, -4*unos, unos, -4*unos],[N-2,N-1,-N+2,-N+1],N,N));

ops.AA = AA;
ops.BB = BB;
ops.CC = CC;
ops.DD = DD;
ops.SS = SS;
ops.AAdec = decomposition(AA,'lu');

if mod(N,2)==0
    k = (2*pi/Lmax) * (-N/2:N/2-1);
else
    k = (2*pi/Lmax) * (-(N-1)/2:(N-1)/2);
end
ops.k = fftshift(k(:));
end
