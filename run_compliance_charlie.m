% Understanding compliance c=Charlie
clear all
clc;
clear all;
clc;
close all;

%% Parameters
period = 1;
n_points=1000;
time_delay=0.1;                         % Delay to get the ventricle start contracting (s)
dt=period/(n_points-1);                     % time step
which_C_lv = "charlie";                 % which capacitance for left ventricle - constant, charlie
Clv_min = 0.006753*1333.22; %ml/mmHg
Clv_max = 0.039*1333.22;    % ml/mmHg
filling_time = 0.8*period;   %  time it takes to fill the LV (s) ??? 
contraction_duration = period - filling_time; % duration of contractile part of the LV
tauS = 0.1*period;          % the smaller, the smaller the derivative in the decay. meaning more curved
tauD = 0.1*period;          % the smaller, the smaller the derivative in the growth curve. meaning more curved
d_C_lv = zeros(1,n_points);
t = linspace(0,period,n_points);

if(which_C_lv=="charlie")
    Clv_function2 = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
    d_C_lv(1) = Clv_function2(mod(t(1)-time_delay,period)); 
elseif(which_C_lv=="constant")
    d_C_lv(1) = Clv_max; 
end

%% Calculate P_mv and P_lv
for(i=2:n_points) %n_points
    % LV capacitance
    if(which_C_lv=="charlie")
        tmod = mod(t(i) -time_delay, period);
        Clv_current = Clv_function2(tmod);
        d_C_lv(i) = Clv_current;
    end
% end
% tmod = mod(t -time_delay, period);
% Clv_current = Clv_function2(tmod);
% d_C_lv = Clv_current;

plot(t, d_C_lv);
xlabel("time(s)");
ylabel("C_{lv}")
title("Capacitance Left Ventricle")

