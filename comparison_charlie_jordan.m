t0=0;                                          % initial simulation time
period = 0.8;                                  % cardiac cycle period (s)
Clv_min = 10;                                   %1.35;%*1333.22; %ml/mmHg
Clv_max = 100;                                  %0.53 * 0.039*1333.22;% 0.65 * 0.039*1333.22;    % ml/mmHg
filling_time = 0.6*period;                     %  time it takes to fill the LV (s) ??? 
contraction_duration = period - filling_time;  % duration of contractile part of the LV
time_delay = 0.01;%0.1;                              % Delay to get the ventricle start contracting (s)
tauS = 0.03*period;                              % the smaller, the smaller the derivative in the decay. meaning more curved
                                                % it also controls how steep P_ao and P_lv the curves going up in systole are 
tauD = 0.05*period;                             % the smaller, the smaller the derivative in the growth curve. meaning more curved
                                              
t = linspace(0,0.8,1000);
d_C_lv = zeros(1,1000);
plot(t, 1./elastance_jordan(t));
hold on;


for i = 1:1000
    tmod = mod(t(i) - time_delay, period);
    d_C_lv(i) = Clv_function_Charlie(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period); %0.005*Clv_max, 0.3*Clv_min
end

plot(t, d_C_lv)
