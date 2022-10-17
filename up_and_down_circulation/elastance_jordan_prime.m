function E_prime = elastance_jordan_prime(t)
    [ncycles, n_points_per_cycle, n_points, dt, which_C_lv, P_pv, d_R_pv, d_C_la, d_R_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao, max_prod] = parameters_circulation;
    tau1 = 0.0725*0.8;
    tau2 = 0.4503*0.8;
    m1 = 2.7463;
    m2 = 21.5683;
    Emin = 0.01;   %mmHg/mL
    Emax = 0.1191; %mmHg/mL

    a=tau1;
    b=tau2;
    c=Emin;
    d=Emax;
    m=m1;
    n=m2;

    prod = (m*(t/a).^(m - 1))./(a*((t/a).^m + 1).*((t./b).^n + 1)) - (m*(t./a).^m.*(t/a).^(m - 1))./(a*((t./a).^m + 1).^2.*((t/b).^n + 1)) - (n*(t/a).^m.*(t/b).^(n - 1))./(b.*((t/a).^m + 1).*((t./b).^n + 1).^2);
    k = (d - c)/max_prod;
    E_prime = k*prod+Emin;
end