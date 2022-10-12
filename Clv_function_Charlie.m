function compliance = Clv_function_Charlie(t, tauS, tauD, TS, CD, CS, T)
    a = tauS;
    b = tauD;
    c = TS;
    d = CD;
    e = CS;
    f = T;
   if (t < c)   
      compliance = d*(e/d).^((1-exp(-t/a))/(1-exp(-c./a)));
   else
      compliance = e*(d/e).^((1-exp(-(t-c)/b))/(1-exp(-(f-c)/b)));
   end
end