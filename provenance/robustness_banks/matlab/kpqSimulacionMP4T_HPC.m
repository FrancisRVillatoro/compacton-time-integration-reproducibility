function result = kpqSimulacionMP4T_HPC(meta,Uinic,outmat)
% M3: fourth-order composition of implicit midpoint.

p = meta.p; q = meta.q;
dx = meta.dx; dt = meta.dt; Lmax = meta.Lmax; Tmax = meta.Tmax;
epsl = meta.epsl; c0 = meta.c0;

x = (0:dx:Lmax)'; N = length(x);
ops = absK_make_operators(N,dx,Lmax);
AA = ops.AA; BB = ops.BB; CC = ops.CC; DD = ops.DD;

Un = Uinic(:);
nmax = round(Tmax/dt);
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

beta1 = (2+2^(1/3)+2^(-1/3))/3;
beta2 = 1 - 2*beta1;
beta3 = beta1;
betas = [beta1 beta2 beta3];

newtonStats.count = 0;
newtonStats.sum = 0;
newtonStats.max = 0;
newtonStats.fail_count = 0;

status = "ok";
t_blowup = NaN;

tic_wall = tic;
cpu0 = cputime;

linearOp = -c0*BB + epsl*DD;

for nt = 1:nmax
    t = nt*dt;

    for s = 1:3
        beta = betas(s);

        ApBD = AA + beta*dt/2*linearOp;
        AmBDtUn = (AA - beta*dt/2*linearOp)*Un;

        TOL = meta.newton_tol;
        Ukp1 = Un;
        Uk = zeros(size(Un));

        newtIter = 1;
        while (norm(Uk-Ukp1) > TOL) && (newtIter < meta.newton_max_iter)
            Uk = Ukp1;
            Ukn = (Uk+Un)/2;

            Dp = spdiags(abs(Ukn).^(p-1),0,N,N);
            Dq = spdiags(abs(Ukn).^(q-1),0,N,N);

            Apq = ApBD + beta*dt/2*BB*p*Dp + beta*dt/2*CC*q*Dq;

            Bpk = beta*dt*BB*(abs(Ukn).^(p-1).*Ukn - p/2*abs(Ukn).^(p-1).*Uk);
            Cqk = beta*dt*CC*(abs(Ukn).^(q-1).*Ukn - q/2*abs(Ukn).^(q-1).*Uk);

            Ukp1 = Apq\(AmBDtUn-Bpk-Cqk);
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
    end

    if status ~= "ok"
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
