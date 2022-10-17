function f = rhs_circuit(t, v)
    % get global parameters - necessary because function should depend on t
    % and one other input variable
   [ncycles, n_points_per_cycle, n_points, dt, which_C_lv, P_pv, R_pv, d_C_la, R_mv, R_mv_closed, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av,R_av_closed, R_ao, d_C_ao, max_prod] = parameters_circulation;
    % Calculate capacitance
    if(which_C_lv=="charlie")
        tmod = mod(t -time_delay, period);
        func_Clv = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
        Clv = func_Clv(tmod); %0.005*Clv_max, 0.3*Clv_min
        func_Clv_prime = @(t) Clv_function_Charlie_prime(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
        Clv_prime = func_Clv_prime(tmod);
    elseif(which_C_lv =="elastance_jordan")
        tmod = mod(t -time_delay, period);
         Clv = 1/elastance_jordan(tmod);
         Clv_prime = -elastance_jordan_prime(tmod)/(elastance_jordan(tmod))^2;
    else
        error("Clv function not implemented")
    end
    % Identify variables so it is easy to debug
    d_P_la = v(1);
    d_P_lv = v(2);
    d_P_ao = v(3);
    P_pv_constant = P_pv(1);
    
    d_Q_pv = ( P_pv_constant - d_P_la)/R_pv;
    % Opening and closing MV
    if(t==t0 || d_P_la > d_P_lv ) % open
        d_R_mv = R_mv;
        %d_Q_mv = (d_P_la - d_P_lv)/d_R_mv;
    else                          % closed
        d_R_mv = R_mv_closed;
        %d_Q_mv = 0;
    end
    d_Q_mv = (d_P_la - d_P_lv)/d_R_mv;

    % Opening and closing AV
    if(t==t0 || d_P_ao >= d_P_lv ) % closed
        d_R_av = R_av_closed;
%         d_Q_av = 0;
    else                          % open
        d_R_av = R_av;
%         d_Q_av = (d_P_lv - d_P_ao) / d_R_av;
    end
    d_Q_av = (d_P_lv - d_P_ao) / d_R_av;
    ratio = -Clv_prime/Clv;
    f(1,1) = 1/d_C_la * (d_Q_pv - d_Q_mv);
    f(2,1) =  ratio * d_P_lv - 1/Clv * d_Q_av + 1/Clv *d_Q_mv;
    f(3,1) = 1/d_C_ao * (- d_P_ao/R_ao + d_Q_av );
end