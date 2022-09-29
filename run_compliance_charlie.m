% Understanding compliance c=Charlie
clear all
clc;

% period = 0.8                 %  period of cardiac cycle (s)
% filling_time = 0.8*period   %  time it takes to fill the LV (s) ??? 
% time_delay = 0.0            %  when LV starts contracting(s)
% contraction_duration = period - filling_time; % duration of contractile part of the LV
% tauS = 0.1*period;          % the smaller, the smaller the derivative in the decay. meaning more curved
% tauD = 0.4*period;          % the smaller, the smaller the derivative in the growth curve. meaning more curved
% Clv_min = 0.006753;         % minimum LV compliance
% Clv_max = 0.039;            % maximum LV compliance
% Clv_function2 = @(t) Clv_function_Charlie(t, tauS, tauD, contraction_duration, Clv_max, Clv_min, period);
% t=linspace(0,0.8,1000);
% dt=0.001;
% d_C_lv=zeros(1,1000);
% 
%       
% for i=1:1000
%     tmod = mod(t(i) - time_delay, period);
%     tmod_future = mod(t(i) + dt - time_delay, period);
%     d_C_lv(1,i) = Clv_function2(tmod); 
% end
% plot(t, d_C_lv);



clear all;
clc;
%close all;

%% Pressure data

which_P = "mmhg_wiggers";   
[t, pp1_MV,pp1_LV,pp1_AO] = interpolation_pressure(which_P);

%% Parameters

period = 1%0.8; % cardiac cycle period (s)
n_points=1000;
av_shut = 1;                            % initial AV diode state
time_delay=0.1;                         % Delay to get the ventricle start contracting (s)
dt=period/n_points;                     % time step
which_C_lv = "charlie";                 % which capacitance for left ventricle - constant, charlie
R_av_reference = 0.04;                  % Resistance open AV - based on Moyer thesis
R_av_closed = 1e+10;                    % Resistance to represent closed diode    
d_R_ao=0.49;
P_ao_initial=81.8;
P_lv_initial=3.67;
d_C_ao = 1.47;
Clv_min = 0.006753*1333.22; %ml/mmHg
Clv_max = 0.039*1333.22;    % ml/mmHg
filling_time = 0.8*period;   %  time it takes to fill the LV (s) ??? 
contraction_duration = period - filling_time; % duration of contractile part of the LV
tauS = 0.1*period;          % the smaller, the smaller the derivative in the decay. meaning more curved
tauD = 0.1*period;          % the smaller, the smaller the derivative in the growth curve. meaning more curved

%% Initialization
subplot(2,3,2)
load("mv_flow_prisco")   
t_MV1=(mv_flow_prisco(54:86,1)-mv_flow_prisco(54,1))*0.1/(mv_flow_prisco(86,1)-mv_flow_prisco(54,1));
q_MV1=mv_flow_prisco(54:86,2);
t_MV3_aux=[((mv_flow_prisco(94:98,1)-0.4783))' (mv_flow_prisco(1:53,1)+(mv_flow_prisco(98,1)-0.5))'];
t_MV3 = t_MV3_aux(1,1) + (t_MV3_aux-t_MV3_aux(1,1))*0.3./(t_MV3_aux(1,58)-t_MV3_aux(1,1));
q_MV3=[(mv_flow_prisco(94:98,2))' (mv_flow_prisco(1:53,2))'];
xq1=linspace(0,0.1,125);
xq2=linspace(0.1,0.5,500);
xq3=linspace(0.57,period,375);
qq1_MV1 = interp1(t_MV1,q_MV1,xq1);
qq1_MV2= 0*xq2;
qq1_MV3 = interp1(t_MV3,q_MV3,xq3);
%t=[xq1 xq2 xq3];
d_Q_mv=[qq1_MV1 qq1_MV2 qq1_MV3];
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
        tmod_past = mod(t(i) - dt - time_delay, period);
        Clv_current = Clv_function2(tmod);
        Clv_past = Clv_function2(tmod_past);
        d_C_lv(i-1) = Clv_past;
        d_C_lv(i) = Clv_current;
    end
end


plot(t, d_C_lv);
xlabel("time(s)");
ylabel("C_{lv}")
title("Capacitance Left Ventricle")

