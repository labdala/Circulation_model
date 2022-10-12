function f = rhs_circuit(t, v)
    % get global parameters - necessary because function should depend on t
    % and one other input variable
   [ncycles, n_points_per_cycle, n_points, dt, which_C_lv, which_Q_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao] = parameters;
    % Calculate capacitance
    if(which_C_lv=="charlie")
        tmod = mod(t -time_delay, period);
        func_Clv = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
        Clv = func_Clv(tmod); %0.005*Clv_max, 0.3*Clv_min
        func_Clv_prime = @(t) Clv_function_Charlie_prime(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
        Clv_prime = func_Clv_prime(tmod);
    else
        error("Clv function not implemented")
    end
    % Identify variables so it is easy to debug
    d_P_lv = v(1);
    d_P_ao = v(2);
    vol_lv = v(3);
   
    Q_mv = mv_inflow(t);
    if(t==t0 || d_P_ao > d_P_lv )
        R_av = R_av_closed;
        d_Q_av = 0;
    else
        R_av = R_av_reference;
        d_Q_av = (d_P_lv - d_P_ao) / R_av_reference;
    end
    ratio = -Clv_prime/Clv;

    f(1,1) =  ratio * d_P_lv - 1/Clv * d_Q_av + 1/Clv *Q_mv;
    f(2,1) = 1/d_C_ao * (- d_P_ao/d_R_ao + d_Q_av );
    f(3,1) = (Q_mv - d_Q_av );
end