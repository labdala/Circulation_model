close all;
clear all;
clc;
load("skagseth2_flow_pv");
load("skagseth2_P_la");
t_la= skagseth2_P_la(:,1);
P_la=skagseth2_P_la(:,2);
t_pv=skagseth2_flow_pv(:,1);
Q_pv=skagseth2_flow_pv(:,2)*1000/60;
t=linspace(0,1,1000);
pp1_la = interp1(t_la,P_la,t);
qq1_pv = interp1(t_pv,Q_pv,t);
qq1_pv(1000)= Q_pv(68);
P_l= 17*ones(1000,1)'; %mmHg

vol_inflow=0; %ml/s

for i=2:1000
    vol_inflow = vol_inflow + (qq1_pv(i-1)+qq1_pv(i))/2*0.001;
end


%3 cardiac cycles
vol_inflow=vol_inflow/3;

% to make about 5L/min or 70ml/s
vol_inflow=vol_inflow*10


subplot(2,1,1)
plot(t, qq1_pv);
hold on;
plot(skagseth2_flow_pv(:,1), skagseth2_flow_pv(:,2)*1000/60, 'o')
xlabel("time")
ylabel("SPV Flow (ml/s)")
xlim([0,1]);

R0=sum((P_l-pp1_la).*qq1_pv,"all");
R0=R0/sum(qq1_pv.*qq1_pv, "all");

subplot(2,1,2)
plot(t, pp1_la);
hold on;
plot(skagseth2_P_la(:,1), skagseth2_P_la(:,2), 'o')
plot(t, P_l-qq1_pv*R0);
xlabel("time")
ylabel("LA Pressure (mmHg)")
xlim([0,1]);
legend("Pressure data", "Data points","Approximated pressure with approx. resistance")

% this LA has 5 pulmonary veins - calculating resistance for each:
R_pv=R0/5

figure
plot(t, (P_l-pp1_la)./qq1_pv)
hold on;
plot(t, R0*ones(1000,1))
plot(t, mean((P_l-pp1_la)./qq1_pv)*ones(1000,1))
title("R (mmHg*s)/ml")
legend("R", "R_0", "R_average")

