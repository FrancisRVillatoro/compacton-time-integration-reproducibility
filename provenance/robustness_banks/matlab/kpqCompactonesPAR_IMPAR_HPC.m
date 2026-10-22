function Uinic = kpqCompactonesPAR_IMPAR_HPC(p,q,dx,Lmax,vel,x0,t)
% Periodic compacton initial condition for |K|(p,p).
% q is kept for compatibility with the old calling convention.

x = (0:dx:Lmax)';
n = p;
a = 1;
b = 1;
c = vel;

alpha = 2*n/(n+1)*(abs(c)/a);
beta  = (n-1)/(2*n)*sqrt(a/b);
gamma = 1/(n-1);

phase = x - (sign(c)^(n-1))*c*t - x0;
width = pi/(2*beta);

Uinic = zeros(size(x));

for m = 0:4
    ph = phase + m*(Lmax+dx);
    mask = abs(ph) <= width;
    Uinic = Uinic + sign(c)*alpha^gamma*cos(ph*beta).^(2*gamma).*mask;
end
end
