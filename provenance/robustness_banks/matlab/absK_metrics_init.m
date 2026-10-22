function metrics = absK_metrics_init(nsave)
metrics.t = nan(nsave,1);
metrics.I2 = nan(nsave,1);
metrics.rho_hw_1_2 = nan(nsave,1);
metrics.rho_hw_2_3 = nan(nsave,1);
metrics.rho_hw_3_4 = nan(nsave,1);
metrics.min_u = nan(nsave,1);
metrics.max_u = nan(nsave,1);
metrics.max_abs_u = nan(nsave,1);
metrics.isave = 0;
end
