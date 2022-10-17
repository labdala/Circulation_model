% test backward_euler_vec
% y''+y=0
clear all;
clc;
close all;

f = @(t,v) [v(2); -v(1)];
ic = [1, 0];
tspan = [0,20];
[t,v]=ode45(f, tspan, ic);
plot(t,v(:,1))
hold on;

nsteps = 120;
[t_be, v_be] = backward_euler_vec(f, tspan, ic , nsteps);
plot(t_be,v_be(:,1))
hold on;
nsteps = 5000;
[t_be, v_be] = backward_euler_vec(f, tspan, ic , nsteps);
plot(t_be,v_be(:,1));
hold on;

nsteps = 120;
[t_fe, v_fe] = forward_euler_vec(f, tspan, ic , nsteps);
plot(t_fe,v_fe(:,1))
hold on;

nsteps = 5000;
[t_fe, v_fe] = backward_euler_vec(f, tspan, ic , nsteps);
plot(t_be,v_be(:,1))

legend('ode45', 'BE, nsteps=120', 'BE, nsteps=5000', 'FE, nsteps=120', 'FE, nsteps=5000')
xlim(tspan)