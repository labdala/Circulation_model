function E = elastance_jordan(t)
    tau1 = 0.0725*0.8;
    tau2 = 0.4503*0.8;
    m1 = 2.7463;
    m2 = 21.5683;
    Emin = 0.01;   %mmHg/mL
    Emax = 0.1191; %mmHg/mL

    g1 = (t./tau1).^m1;
    g2 = (t./tau2).^m2;

    prod = (g1./(1+g1)).*(1./(1+g2));
    k = (Emax - Emin)/max(prod);
    E = k*prod+Emin;
end