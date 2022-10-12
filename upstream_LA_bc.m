clear all;
clc;
close all;
%% Parameters
P_l=17; % 22 %mmHg
which_Q = "data_broome_MVopen_passive";        % flux type - piecewise_quadratic, data_broome, data_broome_MVopen_passive, data_broome_MVopen_active, data_prisco
R_pv = 0.4; % 2 % mmHg · s/mL

[t_la,P_la_pironet, P_0, P_1]=interpolation_pressure("mmhg_pironet_MVopen_passive");
%% Plot input data 
[t,Q,Q_AO,Q_PV]=interpolation_flux(which_Q);
subplot(1,3,1)
plot(t,Q,':.');
hold on;
plot(t,Q_PV,':.');
plot(t,Q_AO,':.');
legend("MV flow","PV inflow", "AO inflow")
title("input data")
xlabel("t(s)")
ylabel("Q (ml/s)")

%% Flux assign and plot
subplot(1,3,2)
plot(t,Q_PV,':.');
title("LA inflow from pulmonary veins: known -  from literature - IB model should give us this info")
xlabel("t(s)")
ylabel("Q_{PV} (ml/s)")

P_la_upstr = zeros(1000,1);
P_la_upstr(1) = P_la_pironet(1);
dt=1/999;
C=-7;
%% Pressure upstream 
for i=1:999
    %P_la_upstr(i+1)= P_l-Q_PV(i)*R_pv;
    P_la_upstr(i+1)= P_la_upstr(i)*(1+dt/C*R_pv)+Q_PV(i)*dt/C;
end



subplot(1,3,3)
plot(t, P_l*ones(1000,1));
hold on;
plot(t, P_la_upstr);
plot(t, P_la_pironet);
legend("P_{lungs}", "P_{la}", "P_{la, literature}");
title("Pressure LA upstream: Boundary condition imposed on PVs")
xlabel("t(s)")
ylabel("P (mmHg)")