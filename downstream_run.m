clear all;
clc;
close all;

test_code = 0;
test_with_comparison = 1;
plot_with_comparison = 0;

%% Parameters
[ncycles, n_points_per_cycle, n_points, dt, which_C_lv, which_Q_mv, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av_reference,R_av_closed, d_R_ao, d_C_ao] = parameters;
%% Initialization
P_ao_initial = 77; %81.8;  % mmHg
P_lv_initial = 8.5; %3.67; % mmHg

tmod0 = mod(0 -time_delay, period);
func_Clv = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
vol_lv_initial = func_Clv(tmod0)*P_lv_initial;
cycle = 1;                                  % initialize the variable that tracks in which cycle we are
tf = period *ncycles;
ic = [P_lv_initial P_ao_initial vol_lv_initial];
nsteps = n_points-1;

[tvec, d_Q_mv] = mv_inflow_vec(which_Q_mv, ncycles, n_points, period, t0);
figure();
%[tp,y] = forward_euler_vec(@rhs_circuit, [t0,tf], ic , nsteps);
[tp,y] = ode45(@rhs_circuit, linspace(t0,tf,n_points), ic , nsteps);
d_P_lv = y(:,1);
d_P_ao = y(:,2);
vol_lv = y(:,3);

%% Tests - give user opportunity to choose to test or not
% % check if number of inductors + number of capacitors  = number of ODEs
% 
% % check if aortic valve is stenotic
% 
%% Post-processing
% Calculate capacitance over time
d_C_lv=zeros(1,length(d_P_ao))';
d_Clv_prime=zeros(1,length(d_P_ao))';
for i = 1:length(d_P_ao)
    tmod = mod(tp(i) - time_delay, period);
    d_C_lv(i) = Clv_function_Charlie(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period); %0.005*Clv_max, 0.3*Clv_min
    d_Clv_prime(i) = Clv_function_Charlie_prime(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
end

% Calculate LV volume over time
% vol_lv = zeros(1,length(tp));
% vol_lv(1) = d_C_lv(1) * d_P_lv(1);
% for  i = 2 : length(tp)
%     vol_lv(i) = vol_lv(i-1) + dt * (d_Q_mv(i-1)-d_Q_av(i-1));
% end

% % Calculate stroke volume 
% stroke_volume = zeros(1,ncycles);
% for i=1:length(d_Q_av)-1
%     stroke_volume(cycle) = stroke_volume(cycle) + dt * (d_Q_av(i+1)+d_Q_av(i))/2;
% end
% 
% stroke_volume
% 
% cardiac_output = (stroke_volume/1000) /(period/60) %L/min

% Calculate Q_av given the P_ao and P_lv
d_Q_av = zeros(length(tp),1);
pdiff = zeros(length(tp),1);

for i = 1 : length(tp)
    if d_P_lv(i) < d_P_ao(i)
        d_Q_av(i) = 0;
        d_R_av(i) = R_av_closed; 
    else
        d_Q_av(i) = (d_P_lv(i) - d_P_ao(i))/R_av_reference;
        d_R_av(i) = R_av_reference; 
        pdiff(i) = d_P_lv(i) - d_P_ao(i);
    end
end

% Calculate Q_ao
d_Q_ao = d_P_ao./d_R_ao;

% Calculate average MV and AO flux during the cycle
average_Q_mv_cycle = zeros(ncycles,1);
average_Q_ao_cycle = zeros(ncycles,1);

for n=1:ncycles
    initial_i = fix((1 + (n-1)*n_points_per_cycle));
    final_i = fix(n*n_points_per_cycle);
    for i =  initial_i: final_i-1
        average_Q_mv_cycle(n) = average_Q_mv_cycle(n) + dt*(d_Q_mv(i+1) + d_Q_mv(i))/2; 
        average_Q_ao_cycle(n) = average_Q_ao_cycle(n) + dt*(d_Q_ao(i+1) + d_Q_ao(i))/2;
    end
end


% Get d_Q_mv over cycle
[t, d_Q_mv] = mv_inflow_vec(which_Q_mv, ncycles, n_points, period, t0);

% Conversion to L/min
d_Q_ao = d_Q_ao * 0.06;
d_Q_av = d_Q_av * 0.06;
d_Q_mv = d_Q_mv * 0.06;
%vol_lv = vol_lv;

%% Plots pressure, LV capacitance, Diode state

if (test_with_comparison)
    load('murgo_flow_ao_ml_sec.mat')
    load('murgo_pressure_Ao.mat')
    load('murgo_pressure_LV.mat')
end

figure()
% Pressure data
subplot(2,3,1)
which_P = "mmhg_wiggers";   
[tt, pp1_MV,pp1_LV,pp1_AO] = interpolation_pressure(which_P);
load("wiggers_la_pressure.csv")
load("wiggers_ao_pressure.csv")
load("wiggers_lv_pressure.csv")
plot(wiggers_la_pressure(:,1),wiggers_la_pressure(:,2), 'ro')
hold on;
plot(tt, pp1_MV, 'r-')
plot(wiggers_ao_pressure(:,1),wiggers_ao_pressure(:,2), 'go')
plot(tt, pp1_AO, 'g-')
plot(wiggers_lv_pressure(:,1),wiggers_lv_pressure(:,2), 'bo')
plot(tt, pp1_LV, 'b-')
ylabel("P(mmHg)")
legend("LA", "LA interp","AO", "AO interp","LV", "LV interp")
title("Wiggers diagram - for reference")

% MITRAL VALVE FLOW - PRISCO DATA
subplot(2,3,2)
plot(t, d_Q_mv);
xlabel("time (s)")
ylabel("flux MV (L/min)")
title("Input data - flux MV")

% PRESSURE PLOTS
subplot(2,3,3)
length(tp)
plot(tp, d_P_lv, '-r');
hold on;
plot(tp, d_P_ao, '-g');
plot(murgo_pressure_Ao(:,1)-0.29, murgo_pressure_Ao(:,2))
plot(murgo_pressure_LV(:,1)-0.29, murgo_pressure_LV(:,2))
hold off;
xlabel("time(s)");
ylabel("Pressure (mmHg)")
legend("P_{lv}", "P_{ao}", "P_{ao, murgo}", "P_{lv,murgo}")
title("Pressure Left Ventricle and Aorta")
xlim([0 ncycles*period])

% CAPACITANCE PLOTS
subplot(2,3,4)
plot(tp, d_C_lv);
hold on;
%plot(tp, d_Clv_prime);
hold on;
plot(tp, d_C_ao* ones(length(tp)));
xlabel("time(s)");
ylabel("C (mL/mmHg)")
title("Capacitance Left Ventricle and Capacitance aorta")
%legend("C_{lv}","C_{lv}'", "C_{ao}")
legend("C_{lv}", "C_{ao}")
xlim([0 ncycles*period])

subplot(2,3,5)
plot(tp, d_R_av, 'o');
ylabel("R_{av} (mmHg*s/mL)");
xlabel("time(s)");
title("Resistance / Diode aortic valve")
xlim([0 ncycles*period])
hold on;

subplot(2,3,6)
hold on;
plot(tp, d_Q_ao);
plot(tp,d_Q_av);
plot(t, d_Q_mv);
plot(murgo_flow_ao_ml_sec(:,1)-0.29,murgo_flow_ao_ml_sec(:,2) * 0.06)
ylabel("Q (L/min)");
xlabel("time(s)");
legend("Q_{ao}", "Q_{av}", "Q_{mv}", "Q_{ao}_{murgo}")
title("Flux out of LV into Aorta")
xlim([0 ncycles*period])

figure()
subplot(2,2,1)
plot(vol_lv,d_P_lv);
ylabel("P(mmHg)");
xlabel("V(mL)");
title("PV loop for the left ventricle")


subplot(2,2,2)
%plot(tp,d_P_lv);
hold on;
plot(tp,vol_lv);
%legend("P_{lv}", "V_{lv}")
legend("V_{lv}")
ylabel("volume (mL)")
xlabel("t(s)");
title("LV volume vs time")
xlim([0 ncycles*period])

subplot(2,2,3)
plot(tp,pdiff);
ylabel("P_{LV} - P_{Ao} (mmHg)");
xlabel("t(s)");
title("Check if AV is stenotic (P_{diff}<10 mmHg)")
xlim([0 ncycles*period])

subplot(2,2,3)
plot(tp,pdiff);
ylabel("P_{LV} - P_{Ao} (mmHg)");
xlabel("t(s)");
title("Check if AV is stenotic (P_{diff}<10 mmHg)")
xlim([0 ncycles*period])

% subplot(2,2,4)
% plot(1:ncycles, average_Q_mv_cycle - average_Q_ao_cycle, 'o');
% ylabel("V_{mv}-V_{av}(mL)");
% xlabel("t(s)");
% title("Conservation of mass test")
% xlim([0 ncycles*period])
subplot(2,2,4)
plot(tp, d_Q_mv);
hold on;
plot(tp, d_Q_ao);
hold on; 
plot((1:ncycles)*period, average_Q_mv_cycle - average_Q_ao_cycle, 'o-');
grid on;
ylabel("Q_{Ao}, Q_{mv} (L) and Delta Q (mL)");
xlabel("t(s)");
legend("Q_{mv}", "Q_{ao}", "average(Q_{mv})_{cycle} - average(Q_{ao})_{cycle}")
title("Conservation of mass test")
xlim([0 ncycles*period])
ylim([-10 10])

% 
% 
% 
% 
% 
% %%
% % figure
% % plot([t t+0.8 t+1.6], [d_Q_mv d_Q_mv d_Q_mv], 'o');