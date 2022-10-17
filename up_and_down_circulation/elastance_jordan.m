function E = elastance_jordan(t)
[ncycles, n_points_per_cycle, n_points, dt, which_C_lv, P_pv, d_R_pv, d_C_la, d_R_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao, max_prod] = parameters_circulation;
    Emin = 0.01;   %mmHg/mL
    Emax = 0.1191; %mmHg/mL

    prod = prod_func(t);
    k = (Emax - Emin)/max_prod;
    E = k*prod+Emin;
end