tp = linspace(0,5,10000);
E = elastance_jordan(tp);
E_prime = elastance_jordan_prime(tp);
figure()
plot (tp, E)
hold on
%plot(tp, E_prime)