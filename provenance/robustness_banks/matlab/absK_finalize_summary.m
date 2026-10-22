function summary = absK_finalize_summary(meta,metrics,newtonStats,status,t_blowup,wall_time,cpu_time)
t = metrics.t(:);
rho = metrics.rho_hw_1_2(:);
I2 = metrics.I2(:);

valid = ~isnan(t);
t = t(valid);
rho = rho(valid);
I2 = I2(valid);

summary = meta;
summary.status = string(status);
summary.t_blowup = t_blowup;
summary.wall_time = wall_time;
summary.cpu_time = cpu_time;

if isempty(t)
    summary.I2_initial = NaN;
    summary.I2_final = NaN;
    summary.I2_drift_abs = NaN;
    summary.I2_drift_rel = NaN;
    summary.rho_hw_final = NaN;
    summary.rho_hw_max = NaN;
    summary.rho_hw_post_auc = NaN;
else
    summary.I2_initial = I2(1);
    summary.I2_final = I2(end);
    summary.I2_drift_abs = abs(I2(end)-I2(1));
    summary.I2_drift_rel = abs(I2(end)-I2(1))/max(abs(I2(1)),eps);
    summary.rho_hw_final = rho(end);
    summary.rho_hw_max = max(rho);

    tc = meta.post_t0;
    idx = t >= tc;
    if sum(idx) >= 2
        summary.rho_hw_post_auc = trapz(t(idx),rho(idx));
    else
        summary.rho_hw_post_auc = NaN;
    end
end

summary.min_u_final = metrics.min_u(find(~isnan(metrics.min_u),1,'last'));
summary.max_u_final = metrics.max_u(find(~isnan(metrics.max_u),1,'last'));
summary.max_abs_u_final = metrics.max_abs_u(find(~isnan(metrics.max_abs_u),1,'last'));

summary.newton_count = newtonStats.count;
summary.newton_sum = newtonStats.sum;
summary.newton_avg = newtonStats.sum/max(newtonStats.count,1);
summary.newton_max = newtonStats.max;
summary.newton_fail_count = newtonStats.fail_count;
end
