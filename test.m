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
x1=x1/100;
x2=x2/100;
x3=x3/100;
y1=y1/100;
y2=y2/100;
y3=y3/100;

t=[x1 x2 x3];
Q=[y1 y2 y3];
dt=1/1000;

plot(t,Q)
hold on

R_mv=1;
P_lv = Q - exp(t);
P_la = Q*R_mv+P_lv
plot(t,P_lv)
hold on;
plot(t, P_la)
legend ("Q", "P_lv", "P_la")