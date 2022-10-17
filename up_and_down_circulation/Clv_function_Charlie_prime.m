function compliance = Clv_function_Charlie_prime(t, tauS, tauD, TS, CD, CS, T)
    a = tauS;
    b = tauD;
    c = TS;
    d = CD;
    e = CS;
    f = T;
     if (t < TS)   
      compliance = -(d.*exp(-t./a).*log(e./d).*(e./d).^((exp(-t./a) - 1)./(exp(-c./a) - 1)))./(a.*(exp(-c./a) - 1));
     else
       compliance = -(e.*exp((c - t)./b).*log(d./e).*(d./e).^((exp((c - t)./b) - 1)./(exp((c - f)./b) - 1)))./(b.*(exp((c - f)./b) - 1));
     end
end
