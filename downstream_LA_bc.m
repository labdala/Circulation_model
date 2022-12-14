clear all;
clc;
close all;

%% Pressure data
subplot(2,3,1)
which_P = "mmhg_wiggers";   
[t, pp1_MV,pp1_LV,pp1_AO] = interpolation_pressure(which_P);
load("wiggers_la_pressure.csv")
load("wiggers_ao_pressure.csv")
load("wiggers_lv_pressure.csv")
plot(wiggers_la_pressure(:,1),wiggers_la_pressure(:,2), 'ro')
hold on;
plot(t, pp1_MV, 'r-')
plot(wiggers_ao_pressure(:,1),wiggers_ao_pressure(:,2), 'go')
plot(t, pp1_AO, 'g-')
plot(wiggers_lv_pressure(:,1),wiggers_lv_pressure(:,2), 'bo')
plot(t, pp1_LV, 'b-')
ylabel("P(mmHg)")
legend("LA", "LA interp","AO", "AO interp","LV", "LV interp")
title("Wiggers diagram")

%% Parameters

period = 0.8; % cardiac cycle period (s)
n_points=1000;
av_shut = 1;                            % initial AV diode state
time_delay=0.1;                         % Delay to get the ventricle start contracting (s)
dt=period/n_points;                     % time step
which_C_lv = "charlie";                 % which capacitance for left ventricle - constant, charlie
R_av_reference = 0.015;%0.04;       0.04%           % Resistance open AV - based on Moyer thesis
R_av_closed = 1e+10;                    % Resistance to represent closed diode    
d_R_ao=0.8; %0.8%
P_ao_initial=116;%81.8;
P_lv_initial=6.6%3.67;
d_C_ao =2.5% 2.8;%1.47;
Clv_min = 7.5 * 0.006753*10;%*1333.22; %ml/mmHg
Clv_max = 0.4 * 0.039*1333.22;% 0.65 * 0.039*1333.22;    % ml/mmHg
filling_time = 0.66*period;   %  time it takes to fill the LV (s) ??? 
contraction_duration = period - filling_time; % duration of contractile part of the LV
tauS = 0.1*period;          % the smaller, the smaller the derivative in the decay. meaning more curved
tauD =  0.445*period;        % the smaller, the smaller the derivative in the growth curve. meaning more curved
stroke_volume=0;
%% Initialization
subplot(2,3,2)
load("mv_flow_prisco")   
t_MV1=(mv_flow_prisco(54:86,1)-mv_flow_prisco(54,1))*0.1/(mv_flow_prisco(86,1)-mv_flow_prisco(54,1));
q_MV1=mv_flow_prisco(54:86,2);
t_MV3_aux=[((mv_flow_prisco(94:98,1)-0.4783))' (mv_flow_prisco(2:53,1)+mv_flow_prisco(98,1)-0.4783)'];
t_MV3 = t_MV3_aux(1,1) + (t_MV3_aux-t_MV3_aux(1,1)).*0.3./(t_MV3_aux(1,57)-t_MV3_aux(1,1));
q_MV3=[(mv_flow_prisco(94:98,2))' (mv_flow_prisco(2:53,2))'];
xq1=linspace(0,0.1,125);
xq2=linspace(0.1,0.5,500);
xq3=linspace(0.5,period,375);
qq1_MV1 = interp1(t_MV1,q_MV1,xq1);
qq1_MV2= 0*xq2;
qq1_MV3 = interp1(t_MV3,q_MV3,xq3);
t=[xq1 xq2 xq3];
d_Q_mv=[qq1_MV1 qq1_MV2 qq1_MV3];
d_Q_av=zeros(size(d_Q_mv));
d_Q_ao=zeros(size(d_Q_mv));

plot(t, d_Q_mv);
xlabel("time (s)")
ylabel("flux MV (ml/s)")
title("Input data - flux MV")

% AV diode
d_R_av = zeros(n_points,1);
if(av_shut==1)
    d_R_av(1) = R_av_closed;
elseif(av_shut==0)
    d_R_av(1) = R_av_reference;
end

% AO and PV pressures 
d_P_ao = zeros(n_points,1);
d_P_ao(1) = P_ao_initial;  % mmHg
d_P_lv = zeros(n_points,1);
d_P_lv(1) = P_lv_initial;  % mmHg

% LV capacitance
d_C_lv=zeros(n_points,1);

if(which_C_lv=="charlie")
    Clv_function2 = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period); %0.005*Clv_max, 0.3*Clv_min
    d_C_lv(1) = Clv_function2(mod(t(1)-time_delay,period)); 
elseif(which_C_lv=="constant")
    d_C_lv(1) = Clv_max; 
end

% Initialize volume
vol_lv = zeros(1,1000);
vol_lv(1) = d_C_lv(1) * d_P_lv(1);

%d_C_lv=10*d_C_lv;

%% Calculate P_mv and P_lv
for(i=2:n_points) %n_points
    % LV capacitance
    if(which_C_lv=="charlie")
        tmod = mod(t(i) -time_delay, period);
        tmod_past = mod(t(i) - dt - time_delay, period);
        Clv_current = Clv_function2(tmod);
        Clv_past = Clv_function2(tmod_past);
        d_C_lv(i-1) = Clv_past;
        d_C_lv(i) = Clv_current;
%     elseif(which_C_lv=="constant")
%         d_C_lv(i-1) = Clv_max;
%         d_C_lv(i) = Clv_max;
    end

   if(d_Q_mv(i)==0) 
       %fprintf("i= %i, time= %f, Closed MV due to flux=0 \n", i, t(i));
   end

   if(d_P_lv(i-1)<d_P_ao(i-1)) 
       fprintf("i= %i, time= %f, P_LV=%f, P_AO=%f, R_av=%f, Closing AV\n", i, t(i), d_P_lv(i-1), d_P_ao(i-1),d_R_av(i));
       av_shut=1;
       d_R_av(i) = R_av_closed;
   else
      av_shut=0;
       d_R_av(i) = R_av_reference;
       fprintf("i= %i, time= %f, P_LV=%f, P_AO=%f, R_av=%f, Opening AV\n", i, t(i), d_P_lv(i-1), d_P_ao(i-1),d_R_av(i));
   end

    % Update Q_av
    if(av_shut==0)
        d_Q_av(i-1) = (d_P_lv(i-1)-d_P_ao(i-1))/d_R_av(i-1);
    elseif(av_shut==1)
        d_Q_av(i-1) = 0;
    end

    d_P_lv(i) = (d_P_lv(i-1) * d_C_lv(i-1) + dt * (d_Q_mv(i)-d_Q_av(i-1)))/d_C_lv(i);
    d_P_ao(i) = d_P_ao(i-1) + ( dt * ( d_Q_av(i-1) -d_P_ao(i-1)/d_R_ao))/ d_C_ao;

    vol_lv(i) = vol_lv(i-1) + dt * (d_Q_mv(i-1)-d_Q_av(i-1));
%     d_Q_ao(i) = d_Q_av(i) + d_C_ao * (d_P_ao(i) - d_P_ao(i-1))/dt;
 end


%% Plots pressure, LV capacitance, Diode state
subplot(2,3,3)
plot(t, d_P_lv, '-r');
hold on;
plot(t, d_P_ao, '-g');
hold off;
xlabel("time(s)");
ylabel("Pressure (mmHg)")
legend("P_{lv}", "P_{ao}")
title("Pressure Left Ventricle and Aorta")

subplot(2,3,4)
plot(t, d_C_lv);
hold on;
plot(t, d_C_ao* ones(size(t)));
xlabel("time(s)");
ylabel("C")
title("Capacitance Left Ventricle and Capacitance aorta")
legend("C_{lv}", "C_{ao}")

subplot(2,3,5)
plot(t, d_R_av, 'o');
ylabel("R_{av} (mmHg*s/mL)");
xlabel("time(s)");
title("Resistance / Diode aortic valve")

subplot(2,3,6)
%plot(t, d_Q_av, 'o');
hold on;
plot(t, d_P_ao./d_R_ao);
plot(t,d_Q_av);
plot(t, d_Q_mv);
ylabel("Q (mL/s)");
xlabel("time(s)");
legend("Q_{ao}", "Q_{av}", "Q_{mv}")
title("Flux out of LV into Aorta")

figure()
subplot(1,2,1)
plot(vol_lv,d_P_lv);
ylabel("P(mmHg)");
xlabel("V(mL)");
title("PV loop for the left ventricle")

subplot(1,2,2)
plot(t,d_P_lv);
hold on;
plot(t,vol_lv);
legend("P_{lv}", "V_{lv}")
xlabel("t(s)");
title("Time vs Pressure")
% figure
% t_wiggers=wiggers_lv_pressure(:,1);
% d_C_lv_wiggers = zeros(1,64);
% for i= 2:64
%     tmod = mod(t_wiggers(i) + time_delay, period);
%     tmod_past = mod(t_wiggers(i) - dt + time_delay, period);
%     Clv_current = Clv_function2(tmod);
%     Clv_past = Clv_function2(tmod_past);
%     d_C_lv_wiggers(i-1) = Clv_past;
%     d_C_lv_wiggers(i) = Clv_current;
% end
%     plot(d_C_lv_wiggers.*wiggers_lv_pressure(:,2),wiggers_lv_pressure(:,2));