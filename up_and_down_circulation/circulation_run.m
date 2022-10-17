clear all;
clc;
close all;

load('colorblind_colormap.mat')

test_code = 0;
test_with_comparison = 1;
plot_with_comparison = 0;

%% Parameters
[ncycles, n_points_per_cycle, n_points, dt, which_C_lv, P_pv, R_pv, d_C_la, R_mv, R_mv_closed, tauS, tauD, contraction_duration, Clv_max, Clv_min, period, time_delay, t0, R_av,R_av_closed, R_ao, d_C_ao, max_prod] = parameters_circulation;
%% Initialization
P_ao_initial = 77; %81.8;  % mmHg
P_lv_initial = 8.0; %3.67; % mmHg
P_la_initial = 8.1;
vol_lv_initial = 125; %mL
vol_la_initial = 23.4; %mL
tf = period *ncycles;
ic = [P_la_initial P_lv_initial P_ao_initial];
nsteps = n_points-1;

% [tvec, d_Q_mv] = mv_inflow_vec(which_Q_mv, ncycles, n_points, period, t0);
%figure();
%[tp,y] = forward_euler_vec(@rhs_circuit, [t0,tf], ic , nsteps);
[tp,y] = backward_euler_vec(@rhs_circuit, linspace(t0,tf,n_points), ic , nsteps);
%[tp,y] = ode45(@rhs_circuit, linspace(t0,tf,n_points), ic , nsteps);
d_P_la = y(:,1);
d_P_lv = y(:,2);
d_P_ao = y(:,3);

%% Tests - give user opportunity to choose to test or not
% % check if number of inductors + number of capacitors  = number of ODEs
% 
% % check if aortic valve is stenotic
% 
%% Post-processing
% Calculate compliance over time
d_C_lv=zeros(1,length(d_P_ao))';
d_Clv_prime=zeros(1,length(d_P_ao))';
for i = 1:length(d_P_ao)
    if(which_C_lv=="charlie")
        tmod = mod(tp(i) - time_delay, period);
        d_C_lv(i) = Clv_function_Charlie(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period); %0.005*Clv_max, 0.3*Clv_min
        d_Clv_prime(i) = Clv_function_Charlie_prime(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
    elseif(which_C_lv =="elastance_jordan")
         tmod = mod(tp(i) - time_delay, period);
         d_C_lv(i) = 1/elastance_jordan(tmod);
         d_Clv_prime(i) = -elastance_jordan_prime(tmod)/(elastance_jordan(tmod))^2;
    end
end

% Calculate volume 
vol_lv = (125-d_C_lv(1) .* d_P_lv(1)) + d_C_lv .* d_P_lv;
vol_la = (23.4 -d_C_la(1) .* d_P_la(1)) + d_C_la .* d_P_la;


% Calculate Q_av given the P_ao and P_lv
d_Q_mv = zeros(length(tp),1);
d_Q_av = zeros(length(tp),1);
pdiff = zeros(length(tp),1);
mv_open = zeros(length(tp),1);
av_open = zeros(length(tp),1);

for i = 1 : length(tp)
    if(d_P_la(i) > d_P_lv(i) ) % open
        %d_Q_mv(i) = (d_P_la(i) - d_P_lv(i))/d_R_mv;
        d_R_mv = R_mv;
        mv_open(i) = 1; 
    else                          % closed
        d_R_mv = R_mv_closed;
        mv_open(i) = 0;
    end
    if d_P_lv(i) < d_P_ao(i)
        d_R_av = R_av_closed;
        %d_Q_av(i) = 0;
        av_open(i) = 0; 
    else
        d_R_av = R_av;
        %d_Q_av(i) = (d_P_lv(i) - d_P_ao(i))/d_R_av;
        av_open(i) = 1; 
        pdiff(i) = d_P_lv(i) - d_P_ao(i);
    end
    d_Q_mv(i) = (d_P_la(i) - d_P_lv(i))/d_R_mv;
    d_Q_av(i) = (d_P_lv(i) - d_P_ao(i))/d_R_av;
end

% Calculate Q_ao
d_Q_ao = d_P_ao./R_ao;
d_Q_pv = (P_pv - d_P_la)/R_pv;

% Calculate average MV and AO flux during the cycle
average_Q_pv_cycle = zeros(ncycles,1);
average_Q_ao_cycle = zeros(ncycles,1);

for n=1:ncycles
    initial_i = fix((1 + (n-1)*n_points_per_cycle));
    final_i = fix(n*n_points_per_cycle);
    for i =  initial_i: final_i-1
        average_Q_pv_cycle(n) = average_Q_pv_cycle(n) + dt*(d_Q_pv(i+1) + d_Q_pv(i))/2; 
        average_Q_ao_cycle(n) = average_Q_ao_cycle(n) + dt*(d_Q_ao(i+1) + d_Q_ao(i))/2;
    end
end

% Calculate stroke volume 
stroke_volume = zeros(1,ncycles);
for cycle=1:ncycles
    for i=1:length(d_Q_av)-1
        stroke_volume(cycle) = stroke_volume(cycle) + dt * (d_Q_av(i+1)+d_Q_av(i))/2;
    end
end

stroke_volume

cardiac_output = (stroke_volume/1000) /(period/60) %L/min

% Conversion to L/min
d_Q_ao = d_Q_ao * 0.06;
d_Q_av = d_Q_av * 0.06;
d_Q_mv = d_Q_mv * 0.06;
d_Q_pv = d_Q_pv * 0.06;
%vol_lv = vol_lv;

%% Plots pressure, LV compliance, Diode state

if (test_with_comparison)
    load('murgo_flow_ao_ml_sec.mat')
    load('murgo_pressure_Ao.mat')
    load('murgo_pressure_LV.mat')
end

figure()

% MITRAL VALVE FLOW 
subplot(2,3,1)
plot(tp, P_pv, 'color', colorblind(1,:));
xlabel("time (s)")
ylabel("Pressure (mmHg)")
title("Input data - pressure source PV ")
xlim([0 ncycles*period])

% COMPLIANCE PLOTS
subplot(2,3,2)
plot(tp, d_C_lv, 'color', colorblind(8,:));
hold on;
%plot(tp, d_Clv_prime);
hold on;
plot(tp, d_C_ao* ones(length(tp)), 'color', colorblind(2,:));
hold on;
plot(tp, d_C_la* ones(length(tp)), 'color', colorblind(5,:));
xlabel("time(s)");
ylabel("C (mL/mmHg)")
title("Compliances")
%legend("C_{lv}","C_{lv}'", "C_{ao}")
legend("C_{lv}", "C_{ao}", "C_{la}")
xlim([0 ncycles*period])

% PRESSURE PLOTS
subplot(2,3,3)
plot(tp, d_P_lv, 'color', colorblind(1,:));
hold on;
plot(tp, d_P_ao, 'color', colorblind(2,:));
hold on;
plot(tp, d_P_la, 'color', colorblind(4,:));
hold on;
%plot(murgo_pressure_Ao(:,1)-0.29, murgo_pressure_Ao(:,2), 'color', colorblind(4,:))
hold on;
%plot(murgo_pressure_LV(:,1)-0.29, murgo_pressure_LV(:,2), 'color', colorblind(5,:))
hold off;
xlabel("time(s)");
ylabel("Pressure (mmHg)")
legend("P_{lv}", "P_{ao}","P_{LA}")%, "P_{ao, murgo}", "P_{lv,murgo}")
title("Output Pressures")
xlim([0 ncycles*period])

% MITRAL VALVE FLOW
subplot(2,3,4)
plot(tp, d_Q_mv, 'color', colorblind(1,:));
hold on;
plot(tp, d_Q_pv, 'color', colorblind(2,:));
xlabel("time (s)")
ylabel("flux (L/min)")
legend("Q_{mv}", "Q_{pv}")
title("Flux pulmonary vein and mitral valve")
xlim([0 ncycles*period])

subplot(2,3,5)
plot(tp, av_open,'o', 'color', colorblind(5,:));
hold on;
plot(tp, mv_open,'o', 'color', colorblind(6,:));
ylabel("state");
xlabel("time(s)");
title("Diode states: 1=open; 0=closed")
xlim([0 ncycles*period])
legend("AV", "MV")
hold on;

subplot(2,3,6)
hold on;
plot(tp, d_Q_ao, 'color', colorblind(1,:));
plot(tp,d_Q_av, 'color', colorblind(2,:));
plot(tp, d_Q_mv, 'color', colorblind(5,:));
%plot(murgo_flow_ao_ml_sec(:,1)-0.29,murgo_flow_ao_ml_sec(:,2) * 0.06, 'color', colorblind(4,:))
ylabel("Q (L/min)");
xlabel("time(s)");
legend("Q_{ao}", "Q_{av}", "Q_{mv}")%, "Q_{ao}_{murgo}")
title("Flux out of LV into Aorta")
xlim([0 ncycles*period])

figure()
subplot(2,2,1)
plot(vol_lv,d_P_lv, 'color', colorblind(6,:));
ylabel("P(mmHg)");
xlabel("V(mL)");
title("PV loop for the left ventricle")


subplot(2,2,2)
%plot(tp,d_P_lv);
plot(tp,vol_lv, 'color', colorblind(1,:));
%legend("P_{lv}", "V_{lv}")
ylabel("V (mL)")
xlabel("t(s)");
title("LV volume vs time")
xlim([0 ncycles*period])

subplot(2,2,3)
plot(vol_la,d_P_la, 'color', colorblind(12,:));
ylabel("P_{LA} (mmHg)");
xlabel("V (mL)");
title("PV loop for the left atrium")

subplot(2,2,4)
plot(tp,vol_la, 'color', colorblind(9,:));
ylabel("V(mL)")
xlabel("t(s)");
title("LA volume vs time")
xlim([0 ncycles*period])


figure()
subplot(2,1,1)
plot(tp,pdiff, 'color', colorblind(10,:));
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
subplot(2,1,2)
plot(tp, d_Q_pv-d_Q_ao, 'color', colorblind(1,:));
hold on;
% plot(tp, d_Q_ao);
hold on; 
%plot((1:ncycles)*period, average_Q_pv_cycle - average_Q_ao_cycle, 'o-');
grid on;
ylabel("Q_{Ao} -  Q_{pv} (L)");
xlabel("t(s)");
%legend("Q_{mv}", "Q_{ao}", "average(Q_{mv})_{cycle} - average(Q_{ao})_{cycle}")
title("Conservation of mass test (average should be 0)")
xlim([0 ncycles*period])

subtraction_average=average_Q_pv_cycle - average_Q_ao_cycle


% 
% 
% 
% 
% 
% %%
% % figure
% % plot([t t+0.8 t+1.6], [d_Q_mv d_Q_mv d_Q_mv], 'o');