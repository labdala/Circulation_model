d_C_lv_jordan=zeros(1,length(d_P_ao))';
d_C_lv_charlie=zeros(1,length(d_P_ao))';
d_Clv_prime=zeros(1,length(d_P_ao))';
for i = 1:length(d_P_ao)
        tmod = mod(tp(i) - time_delay, period);
        d_C_lv_charlie(i) = Clv_function_Charlie(tmod, tauS, tauD, contraction_duration, Clv_max, Clv_min, period); %0.005*Clv_max, 0.3*Clv_min
        d_C_lv_jordan(i) = 1/elastance_jordan(tmod);
end

plot(tp, d_C_lv_charlie);
hold on;
plot(tp, d_C_lv_jordan);
