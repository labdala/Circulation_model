[ncycles, n_points_per_cycle, n_points, dt, which_C_lv, P_pv, R_pv, d_C_la, R_mv, R_mv_closed, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av,R_av_closed, R_ao, d_C_ao, max_prod] = parameters_circulation;
t= linspace(0,3.2);

%tmod = mod(t -0.1, 0.8);
func_Clv_prime = @(t,c) Clv_function_Charlie_prime(mod(t -0.1, 0.8), tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
        
[t_c,c_curve] = ode45(func_Clv_prime, [0, 3.2], d_C_lv(1));

plot(tp, d_C_lv, 'r');
hold on
plot(t_c, c_curve, 'b');
legend('C_lv', 'integral')

plot(t,func_Clv_prime(t))