function Q_mv = mv_inflow(t)
% TODO: make this efficient
[ncycles, n_points_per_cycle, n_points, dt, which_C_lv, which_Q_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao] = parameters;
[tvec, d_Q_mv] = mv_inflow_vec(which_Q_mv, ncycles, n_points, period, t0);

i = fix(t/dt +1);
Q_mv = d_Q_mv(i);

end