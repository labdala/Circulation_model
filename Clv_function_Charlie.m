function compliance = Clv_function_Charlie(t, tauS, tauD, TS, CD, CS, T);
   if (t < TS)   
      compliance = CD*(CS/CD)^((1-exp(-t/tauS))/(1-exp(-TS/tauS)));
   else
      compliance = CS*(CD/CS)^((1-exp(-(t-TS)/tauD))/(1-exp(-(T-TS)/tauD)));
   end
end