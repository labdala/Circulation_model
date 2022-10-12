function [ncycles, n_points_per_cycle, n_points, dt, which_C_lv, which_Q_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao] = parameters;
    period = 0.8;
    ncycles = 4;
    n_points_per_cycle = 1000;
    n_points = ncycles * n_points_per_cycle;
    dt = period * ncycles/(n_points-1);                         % time step
    which_C_lv = "charlie";                       % which capacitance for left ventricle - constant, charlie
    which_Q_mv = "constant";
    %
    t0=0;                                          % initial simulation time
    period = 0.8;                                  % cardiac cycle period (s)
    Clv_min = 5;% Jordan = 10; %5.5;                                   %1.35;%*1333.22; %ml/mmHg
    Clv_max = 67;% Jordan = 100; %62;                                  %0.53 * 0.039*1333.22;% 0.65 * 0.039*1333.22;    % ml/mmHg
    filling_time = 0.6*period;                     %  time it takes to fill the LV (s) ??? 
    contraction_duration = period - filling_time;  % duration of contractile part of the LV
    time_delay = 0.1;                              % Delay to get the ventricle start contracting (s)
    tauS = 0.03*period; %0.05*period;                            % the smaller, the smaller the derivative in the decay. meaning more curved
                                                   % it also controls how steep P_ao and P_lv the curves going up in systole are 
    tauD = 0.05*period; %0.09*period;                             % the smaller, the smaller the derivative in the growth curve. meaning more curved
                                                   % it also controls how steep P_ao and P_lv the curves going down in systole are 
    
    R_av_reference = 0.007;%Jordan = 0.0043;                         % Resistance open AV - based on Moyer thesis should be 0.04
    R_av_closed = 1e+10;                           % Resistance to represent closed diode    
    d_R_ao=0.9; 
    d_C_ao = 1.47; 
end