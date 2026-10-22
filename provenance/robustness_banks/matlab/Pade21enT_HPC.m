function result = Pade21enT_HPC(meta,Uinic,outmat)
% M2: Padé(2,1) in time, L-stable. Artificial viscosity only.

p = meta.p; q = meta.q;
dx = meta.dx; dt = meta.dt; Lmax = meta.Lmax; Tmax = meta.Tmax;
epsl = meta.epsl; c0 = meta.c0;

x = (0:dx:Lmax)'; N = length(x);
ops = absK_make_operators(N,dx,Lmax);
AA = ops.AA; BB = ops.BB; CC = ops.CC; DD = ops.DD; AAdec = ops.AAdec;

Un = Uinic(:);
nmax = ceil(Tmax/dt);
nplt = max(1,floor((Tmax/meta.nsave_target)/dt));
nsave = floor(nmax/nplt)+2;

if meta.save_snapshots
    Udata = zeros(N,nsave);
    tdata = zeros(1,nsave);
else
    Udata = [];
    tdata = [];
end

metrics = absK_metrics_init(nsave);
metrics = absK_metrics_add(metrics,Un,0,p,dx,ops.k);

isnap = 1;
if meta.save_snapshots
    Udata(:,isnap) = Un;
    tdata(isnap) = 0;
end

newtonStats.count = 0;
newtonStats.sum = 0;
newtonStats.max = 0;
newtonStats.fail_count = 0;

status = "ok";
t_blowup = NaN;

BBB = -c0*BB + epsl*DD;

parcial1 = -AA - dt/4*BBB;
parcial2BB = -dt/4*BB;
parcial2CC = -dt/4*CC;
parcial3 = 2*dt/3*BBB;
parcial4BB = 2*dt/3*BB;
parcial4CC = 2*dt/3*CC;

DepK = zeros(2*N,1);
UVn = zeros(2*N,1);

Abloc = spalloc(2*N,2*N,2*N*(2*5+1));
Abloc(1:2:2*N,2:2:2*N) = AA;
Abloc(2:2:2*N,1:2:2*N) = AA;

Vn = Un + dt/4*(AAdec\(BBB*Un + BB*(abs(Un).^(p-1).*Un) + CC*(abs(Un).^(q-1).*Un)));

UVn(1:2:2*N) = Un;
UVn(2:2:2*N) = Vn;

tic_wall = tic;
cpu0 = cputime;

for nt = 1:nmax
    t = nt*dt;

    AWn = AA*Un - dt/3*(BBB*Un + BB*(abs(Un).^(p-1).*Un) + CC*(abs(Un).^(q-1).*Un));

    TOL = meta.newton_tol;

    Ukp1 = Un;
    Vkp1 = Vn;
    UVkp1 = UVn;
    UVk = zeros(2*N,1);

    newtIter = 1;
    while (norm(UVk-UVkp1) > TOL) && (newtIter < meta.newton_max_iter)
        UVk = UVkp1;
        Uk = Ukp1;
        Vk = Vkp1;

        Dpu = spdiags(abs(Uk).^(p-1),0,N,N);
        Dqu = spdiags(abs(Uk).^(q-1),0,N,N);
        Dpv = spdiags(abs(Vk).^(p-1),0,N,N);
        Dqv = spdiags(abs(Vk).^(q-1),0,N,N);

        Abloc(1:2:2*N,1:2:2*N) = parcial1 + p*parcial2BB*Dpu + q*parcial2CC*Dqu;
        Abloc(2:2:2*N,2:2:2*N) = parcial3 + p*parcial4BB*Dpv + q*parcial4CC*Dqv;

        DepK(1:2:2*N) = (p-1)*parcial2BB*(abs(Uk).^(p-1).*Uk) + ...
                         (q-1)*parcial2CC*(abs(Uk).^(q-1).*Uk);

        DepK(2:2:2*N) = AWn + ...
                         (p-1)*parcial4BB*(abs(Vk).^(p-1).*Vk) + ...
                         (q-1)*parcial4CC*(abs(Vk).^(q-1).*Vk);

        UVkp1 = Abloc\DepK;

        Ukp1 = UVkp1(1:2:2*N);
        Vkp1 = UVkp1(2:2:2*N);

        newtIter = newtIter + 1;
    end

    nit = newtIter - 1;
    newtonStats.count = newtonStats.count + 1;
    newtonStats.sum = newtonStats.sum + nit;
    newtonStats.max = max(newtonStats.max,nit);
    if newtIter >= meta.newton_max_iter
        newtonStats.fail_count = newtonStats.fail_count + 1;
    end

    Un = Ukp1;
    Vn = Vkp1;
    UVn = UVkp1;

    if any(isnan(Un)) || any(isinf(Un))
        status = "nan_or_inf";
        t_blowup = t;
        break
    end

    if max(abs(Un)) > meta.blowup_factor*meta.max_abs_u0
        status = "amplitude_blowup";
        t_blowup = t;
        break
    end

    if mod(nt,nplt)==0 || nt==nmax
        metrics = absK_metrics_add(metrics,Un,t,p,dx,ops.k);
        if meta.save_snapshots
            isnap = isnap + 1;
            Udata(:,isnap) = Un;
            tdata(isnap) = t;
        end
    end
end

wall_time = toc(tic_wall);
cpu_time = cputime - cpu0;

metrics = absK_metrics_trim(metrics);

if meta.save_snapshots
    Udata = Udata(:,1:isnap);
    tdata = tdata(1:isnap);
end

summary = absK_finalize_summary(meta,metrics,newtonStats,status,t_blowup,wall_time,cpu_time);

if meta.save_raw
    save(outmat,'meta','summary','metrics','x','Udata','tdata','-v7.3');
end

result.summary = summary;
result.metrics = metrics;
result.Ufinal = Un;
end
