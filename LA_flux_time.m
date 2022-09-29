%% Flux in the left atrium 
clear all;
close all;
clc;

% x1 = linspace(0,20,200);
% y1 = -0.0145833*x1.^2+1.20833*x1+20.0000;
% x2 = linspace(60,70,100);
% y2 = 7/50*x2.^2-37/2*x2+646;
% x3 = linspace(70,100,300);
% y3 = -13/300*x3.^2+34/5*x3-680/3;
% x1=x1/100;
% x2=x2/100;
% x3=x3/100;
% y1=y1/100;
% y2=y2/100;
% y3=y3/100;

t1 = linspace(0,20,200);
t2 = linspace(20,45,250);
t3 = linspace(45,50,50);
t4 = linspace(50,55,50);
t5 = linspace(55,60,50);
t6 = linspace(60,100,400);
tvec=[t1 t2 t3 t4 t5 t6];
y1 = -0.0145833*x1.^2+1.20833*x1+20.0000;
t2 = linspace(60,70,100);
y2 = 7/50*x2.^2-37/2*x2+646;
t3 = linspace(70,100,300);
y3 = -13/300*x3.^2+34/5*x3-680/3;
t1=x1/100;
t2=x2/100;
t3=x3/100;
y1=y1/100;
y2=y2/100;
y3=y3/100;


t=[x1 x2 x3];
Q=[y1 y2 y3];
dt=1/1000;

figure;
plot(t,Q)
xlabel("time");
ylabel("Q")

figure;
%% Parameters
d_R_mv = zeros(1000,1);
d_R_mv(1) =1e+10; %start closed
d_R_ao = 0.5;
d_C_lv =5; 
%% Initialization
d_P_mv = zeros(1000,1);
d_P_lv = zeros(1000,1);
d_P_mv(1) = 4;  % mmHg
d_P_lv(1) = 14; % mmHg

%% Calculate P_mv and P_lv
for(i=2:1000)
    recalculate_pressures = false;
    d_R_mv(i) = d_R_mv(i-1);
    d_P_lv(i) = d_P_lv(i-1) * d_C_lv * d_R_ao * d_R_mv(i) +  Q(i) * d_R_mv(i) * dt * d_R_ao;
    d_P_mv(i) = d_P_lv(i) + Q(i) * d_R_mv(i) * ( dt * d_R_mv(i) + d_C_lv * d_R_ao * d_R_mv(i)) ;
    
    d_P_lv(i) = d_P_lv(i)/(d_C_lv * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
    d_P_mv(i) = d_P_mv(i)/(d_C_lv * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
    
    % check if the diode stayed in the same mode as previous time step
    % if not, change its state and recalculate the pressures
    %[d_R_mv(i) d_P_mv d_P_lv]

    if(d_R_mv(i)>1 && d_P_lv(i)< d_P_mv(i) )
        fprintf("Open mitral valve\n")
        d_R_mv(i)=1e-5;
        recalculate_pressures = true;
    elseif(d_R_mv(i)<1  && d_P_lv(i) >= d_P_mv(i)) % why is it not entering here?
        fprintf("Close mitral valve\n")
        d_R_mv(i)=1e+10;
        recalculate_pressures= true;
    end
    if(recalculate_pressures)
        fprintf("Recalculating pressures\n")
        d_P_lv(i) = d_P_lv(i-1) * d_C_lv * d_R_ao * d_R_mv(i) +  Q(i) * d_R_mv(i) * dt * d_R_ao;
        d_P_mv(i) = d_P_lv(i) + Q(i) * d_R_mv(i) * ( dt * d_R_mv(i) + d_C_lv * d_R_ao * d_R_mv(i)) ;
        
        d_P_lv(i) = d_P_lv(i)/(d_C_lv * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
        d_P_mv(i) = d_P_mv(i)/(d_C_lv * d_R_ao * d_R_mv(i) + dt * d_R_mv(i));
    end
end



plot(t, d_P_mv, 'ob');
hold on
plot(t, d_P_lv, 'xr');
hold on;
hold off;
xlabel("time(s)");
legend("P_{mv}", "P_{lv}")


figure
plot(t, d_R_mv, 'o');
ylabel("R_{mv}");
xlabel("time(s)");
