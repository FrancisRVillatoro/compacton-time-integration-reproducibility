function metrics = absK_metrics_add(metrics,U,t,p,dx,k)
metrics.isave = metrics.isave + 1;
j = metrics.isave;

Uhat = fft(U);
Etotal = sum(abs(Uhat).^2);
kmax = max(abs(k));

if Etotal == 0
    rho12 = 0;
    rho23 = 0;
    rho34 = 0;
else
    rho12 = sum(abs(Uhat(abs(k) >= 0.5*kmax)).^2)/Etotal;
    rho23 = sum(abs(Uhat(abs(k) >= (2/3)*kmax)).^2)/Etotal;
    rho34 = sum(abs(Uhat(abs(k) >= 0.75*kmax)).^2)/Etotal;
end

metrics.t(j) = t;
metrics.I2(j) = dx*sum(abs(U).^(p+1));
metrics.rho_hw_1_2(j) = rho12;
metrics.rho_hw_2_3(j) = rho23;
metrics.rho_hw_3_4(j) = rho34;
metrics.min_u(j) = min(U);
metrics.max_u(j) = max(U);
metrics.max_abs_u(j) = max(abs(U));
end
