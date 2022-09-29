%% Flux in the left atrium 
clear all;
close all;
clc;

x1 = linspace(0,60,600);
y1 = -0.0145833*x1.^2+1.20833*x1+20.0000;
x2 = linspace(60,70,100);
y2 =  7/50*x2.^2-37/2*x2+646;
x3 = linspace(70,100,300);
y3 = -13/300*x3.^2+34/5*x3-680/3;
x1 = x1/100;
x2 = x2/100;
x3 = x3/100;
y1 = y1/100;
y2 = y2/100;
y3 = y3/100;

t    = [x1 x2 x3];
Q_mv = [y1 y2 y3];
dt   = 1/1000;
R_ao = 1e-5;

%% Parameters
d_R_mv = 0.5;
d_R_ao = 1e+10; %start open
d_C_lv = 5; 
%% Initialization
d_P_mv = 4;  % mmHg
d_P_lv = 14; % mmHg
%f = @(t,x) -0.0145833*t.^2+1.20833*t+20.0000; - x/R_ao;

recalculate_pressures = false;
%% Calculate P_mv and P_lv
tspan  = [0 60];
[t, x] = ode45((@(t,x) -0.0145833*t.^2+1.20833*t+20.0000 - x./R_ao), tspan, 14);
plot(t,x, 'o')
