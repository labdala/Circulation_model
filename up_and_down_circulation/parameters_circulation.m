function [ncycles, n_points_per_cycle, n_points, dt, which_C_lv, P_pv, R_pv, d_C_la, R_mv, R_mv_closed, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av,R_av_closed, R_ao, d_C_ao, max_prod] = parameters_circulation;
    t0=0;                                          % initial simulation time
    period = 0.8;                                  % cardiac cycle period (s)
    ncycles = 4;
    n_points_per_cycle = 1000;
    n_points = ncycles * n_points_per_cycle;
    dt = period * ncycles/(n_points-1);                         % time step
    which_C_lv = "charlie";                       % which capacitance for left ventricle - constant, charlie
    max_prod = max(prod_func(linspace(t0,period*ncycles,n_points)));
    d_C_la = 2;
    P_pv = 17*ones(n_points,1);
    R_pv = 0.6;
    R_mv = 0.004;
    R_mv_closed = 1e+10;  
    %
    Clv_min =0.5;% should be greater than 0.45% Jordan = 10; %5.5;                                   %1.35;%*1333.22; %ml/mmHg
    Clv_max = (75 + 110*Clv_min)/8;%14;% Jordan = 100; %62;                                  %0.53 * 0.039*1333.22;% 0.65 * 0.039*1333.22;    % ml/mmHg
    filling_time = 0.6*period;                     %  time it takes to fill the LV (s) ??? 
    contraction_duration = period - filling_time;  % duration of contractile part of the LV
    time_delay = 0.1;                              % Delay to get the ventricle start contracting (s)
    tauS = 0.03*period; %0.05*period;                            % the smaller, the smaller the derivative in the decay. meaning more curved
                                                   % it also controls how steep P_ao and P_lv the curves going up in systole are 
    tauD = 0.05*period; %0.09*period;                             % the smaller, the smaller the derivative in the growth curve. meaning more curved
                                                   % it also controls how steep P_ao and P_lv the curves going down in systole are 
    
    R_av = 0.01;%Jordan = 0.0043;                         % Resistance open AV - based on Moyer thesis should be 0.04
    R_av_closed = 1e+10;                           % Resistance to represent closed diode    
    R_ao=5; %1.7
    d_C_ao = 2.0; 
end