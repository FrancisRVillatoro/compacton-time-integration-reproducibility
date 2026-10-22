function run_absK_case(cases_csv, task_id)
% Execute one case from a CSV table.

addpath(fileparts(mfilename('fullpath')));

T = readtable(cases_csv, 'TextType','string');
row = T(task_id,:);

bank_dir = char(row.bank_dir);
summary_dir = fullfile(bank_dir,'summary');
raw_dir = fullfile(bank_dir,'raw');

if ~exist(summary_dir,'dir'), mkdir(summary_dir); end
if ~exist(raw_dir,'dir'), mkdir(raw_dir); end

meta = struct();

meta.run_id = row.run_id;
meta.bank = row.bank;
meta.method = row.method;
meta.scenario = row.scenario;

meta.p = row.p;
meta.q = row.p;
meta.a = 1;
meta.b = 1;

meta.dx = row.dx;
meta.dt = row.dt;
meta.Lmax = row.Lmax;
meta.Tmax = row.Tmax;
meta.c0 = row.c0;
meta.c1 = row.c1;
meta.c2 = row.c2;
meta.x01 = row.x01;
meta.x02 = row.x02;

meta.mu = row.mu;
meta.epsl = meta.mu * meta.dx^2;

meta.post_t0 = row.post_t0;

meta.save_raw = logical(row.save_raw);
meta.save_snapshots = logical(row.save_snapshots);
meta.nsave_target = row.nsave_target;

meta.newton_tol = row.newton_tol;
meta.newton_max_iter = row.newton_max_iter;
meta.blowup_factor = row.blowup_factor;

U1 = kpqCompactonesPAR_IMPAR_HPC(meta.p,meta.q,meta.dx,meta.Lmax,meta.c1,meta.x01,0);
U2 = kpqCompactonesPAR_IMPAR_HPC(meta.p,meta.q,meta.dx,meta.Lmax,meta.c2,meta.x02,0);
Uinic = U1 + U2;

meta.max_abs_u0 = max(abs(Uinic));

run_name = sprintf('run_%06d_%s_p%.12g_mu%.0e_dt%.12g_%s', ...
    meta.run_id, char(meta.method), meta.p, meta.mu, meta.dt, char(meta.scenario));

outmat = fullfile(raw_dir, [run_name '.mat']);
outcsv = fullfile(summary_dir, [run_name '.csv']);

fprintf('Running %s\n', run_name);
fprintf('epsl = %.16e = mu*dx^2, mu = %.16e, dx = %.16e\n', meta.epsl, meta.mu, meta.dx);

switch char(meta.method)
    case 'M1'
        result = deFrutosAbsDissipative_HPC(meta,Uinic,outmat);
    case 'M2'
        result = Pade21enT_HPC(meta,Uinic,outmat);
    case 'M3'
        result = kpqSimulacionMP4T_HPC(meta,Uinic,outmat);
    otherwise
        error('Unknown method: %s', char(meta.method));
end

absK_write_summary_csv(result.summary,outcsv);

fprintf('Finished run_id=%d, status=%s\n', meta.run_id, result.summary.status);
end
