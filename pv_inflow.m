function Q_la = pv_inflow(t)
% TODO: make this efficient
[ncycles, n_points_per_cycle, n_points, dt, which_C_lv, which_P_pv, d_R_pv, d_C_la, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao] = parameters_circulation;
[tvec, d_Q_la] = pv_inflow_vec(which_Q_pv, ncycles, n_points, period, t0);

i = fix(t/dt +1);
Q_la = d_Q_la(i);

end