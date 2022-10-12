function d_Q_av = av_flux_valve(P_lv, P_ao)
     [ncycles, n_points, dt, which_C_lv, which_Q_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao] = parameters;
    if(P_lv < P_ao)
        d_Q_av = 0;
    else
        d_Q_av = (P_lv - P_ao) / R_av_reference;
    end
end