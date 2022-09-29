%% Flux in the left atrium 
clear all;
close all;
clc;

%% Parameters
plot_flux_separately = true;            % plot figure with fluxes
which_Q = "data_broome_MVopen_passive";                % flux type - piecewise_quadratic, data_broome, data_broome_MVopen_passive, data_broome_MVopen_active, data_prisco
mv_shut = 0;                            % initial MV diode state
time_delay=0%0.1;                          % NOT SURE WHAT THIS IS
dt=0.001;                               % time step
which_C_lv = "charlie";                 % which capacitance for left ventricle - constant, charlie
R_mv_reference = 0.02;                  % Resistance open MV - based on Moyer thesis
R_mv_closed = 1e+10;                    % Resistance to represent closed diode    

%% Flux assign and plot
[t,Q,Q_AO,Q_PV]=interpolation_flux(which_Q);

if(plot_flux_separately)
    figure
    subplot(3,1,1)
    plot(t,Q,':.');
    title("MV inflow")
    xlabel("t(s)")
    ylabel("Q (ml/s)")
    
    subplot(3,1,2)
    plot(t,Q_PV,':.');
    title("PV inflow")
    xlabel("t(s)")
    ylabel("Q (ml/s)")
    
    subplot(3,1,3)
    plot(t,Q_AO,':.');
    title("AO inflow")
    xlabel("t(s)")
    ylabel("Q (ml/s)")
end 

figure
subplot(2,3,1)
plot(t,Q,':.');
hold on;
plot(t,Q_PV,':.');
plot(t,Q_AO,':.');
title("Prescribed Flux - given by IB model")
xlabel("t(s)")
ylabel("Q (ml/s)")
legend("Flux Mitral Valve","Flux Aortic Valve", "Flux Pulmonary")

%% Initialization

% MV diode
d_R_mv = zeros(1000,1);
if(mv_shut==1)
    d_R_mv(1) = R_mv_closed;
elseif(mv_shut==0)
    d_R_mv(1) = R_mv_reference;
end

% MV pressure 
d_P_mv = zeros(1000,1);
d_P_mv(1) = 4;  % mmHg


% LV capacitance
d_C_lv=zeros(1000,1);
Clv_min = 0.006753;
Clv_max = 0.039;    
if(which_C_lv=="charlie")
    period = 1;                  %  period of cardiac cycle (s)
    filling_time = 0.8*period;   %  time it takes to fill the LV (s) ??? 
    contraction_duration = period - filling_time; % duration of contractile part of the LV
    tauS = 0.1*period;          % the smaller, the smaller the derivative in the decay. meaning more curved
    tauD = 0.4*period;          % the smaller, the smaller the derivative in the growth curve. meaning more curved
    Clv_function2 = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
    d_C_lv(1) = Clv_function2(t(1)+time_delay); 
elseif(which_C_lv=="constant")
    d_C_lv(1) = Clv_max; 
end

% LV pressure
d_P_lv = zeros(1000,1);
d_P_lv(1) = 14; % mmHg


% AO
d_R_ao = 0.5;

%% Calculate P_mv and P_lv
for(i=2:1000) %1000

    % LV capacitance
    if(which_C_lv=="charlie")
        tmod = mod(t(i) + time_delay, period);
        tmod_past = mod(t(i) - dt + time_delay, period);
        Clv_current = Clv_function2(tmod);
        Clv_past = Clv_function2(tmod_past);
        d_C_lv(i-1) = Clv_past;
        d_C_lv(i) = Clv_current;
    elseif(which_C_lv=="constant")
        d_C_lv(i-1) = Clv_max;
        d_C_lv(i) = Clv_max;
    end

    % Initialize flag
    recalculate_pressures = false;
    
    % Assume current mitral valve state remains the same as last time step
    d_R_mv(i) = d_R_mv(i-1);

    % Update MV and LV pressures
    if(mv_shut==0)
        d_P_lv(i) = d_P_lv(i-1) * d_C_lv(i-1) * d_R_ao * d_R_mv(i) +  Q(i) * d_R_mv(i) * dt * d_R_ao;
        d_P_mv(i) = d_P_lv(i) + Q(i) * d_R_mv(i) * ( dt * d_R_mv(i) + d_C_lv(i) * d_R_ao * d_R_mv(i)) ;
        
        d_P_lv(i) = d_P_lv(i)/(d_C_lv(i) * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
        d_P_mv(i) = d_P_mv(i)/(d_C_lv(i) * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
    elseif(mv_shut==1)
        %d_P_lv(i) = 25;
        %d_P_mv(i) = 25;  
        
        d_P_lv(i) =d_P_lv(i-1)*(d_C_lv(i-1)-dt/d_R_ao)/ d_C_lv(i);
        d_P_mv(i) = d_P_lv(i); 

    end

    %[d_R_mv(i) d_P_mv d_P_lv]
    fprintf("i=%i p_mv(i)=%f \n", i, d_P_mv(i));
    
    % check diode state, change it if necessary and recalculate the
    % pressures - open based on flux and close based on pressure???????
    if(mv_shut==1 && Q(i)>0 )
        fprintf("Open mitral valve for i=%i at time t=%f \n", i, t(i))
        mv_shut=0;
        d_R_mv(i)=R_mv_reference;
        recalculate_pressures = true;
    elseif (mv_shut==0  && Q(i)<=0)
        fprintf("Close mitral valve for i=%i at time t=%f \n", i, t(i))
        mv_shut=1;
        d_R_mv(i)=R_mv_closed;
        recalculate_pressures= true;
    end
    
    % recalculate pressures if diode changed state
    if(recalculate_pressures && mv_shut==0) % MV open
        fprintf("Recalculating pressures for open valve\n")
        fprintf("INSIDE RECALC - i=%i,  p_mv(i)=%f \n", i, d_P_mv(i));
        d_P_lv(i) = d_P_lv(i-1) * d_C_lv(i) * d_R_ao * d_R_mv(i) +  Q(i) * d_R_mv(i) * dt * d_R_ao;
        d_P_mv(i) = d_P_lv(i) + Q(i) * d_R_mv(i) * ( dt * d_R_mv(i) + d_C_lv(i) * d_R_ao * d_R_mv(i)) ;
        
        d_P_lv(i) = d_P_lv(i)/(d_C_lv(i) * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
        d_P_mv(i) = d_P_mv(i)/(d_C_lv(i) * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
        fprintf("INSIDE RECALC  FIN- i=%i,  p_mv(i)=%f \n", i, d_P_mv(i));
    elseif(recalculate_pressures && mv_shut==1) % MV closed
        fprintf("Recalculating pressures for open valve\n")
        fprintf("INSIDE RECALC - i=%i,  p_mv(i)=%f \n", i, d_P_mv(i));
        %d_P_lv(i) = 25;
        %d_P_mv(i) = 25;  
        d_P_lv(i) =d_P_lv(i-1)*(d_C_lv(i-1)-dt/d_R_ao)/ d_C_lv(i);
        d_P_mv(i) = d_P_lv(i); 
    end
end


%% Plots pressure, LV capacitance, Diode state
subplot(2,3,2)
plot(t, d_P_mv, 'ob');
hold on
plot(t, d_P_lv, 'xr');
hold on;
%plot(t, Q, '-');
hold off;
xlabel("time(s)");
ylabel("Pressure (mmHg)")
legend("P_{mv}", "P_{lv}")
title("Pressure Mitral Valve and Left Ventricle")

subplot(2,3,3)
plot(t, d_C_lv);
xlabel("time(s)");
ylabel("C_{lv}")
title("Capacitance Left Ventricle")

subplot(2,3,4)
plot(t, d_R_mv, 'o');
ylabel("R_{mv} (mmHg*s/mL)");
xlabel("time(s)");
title("Resistance / Diode mitral valve")

subplot(2,3,5)
plot(t, d_P_lv./d_R_ao, 'o');
ylabel("Q_{lv} (mL/s)");
xlabel("time(s)");
title("Flux out of LV into Aorta")
