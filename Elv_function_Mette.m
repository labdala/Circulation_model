function elastance = Elv_function_Mette(t, tce, th, Emin, Emax)
aphi = 0.9;
bphi = 0.25;
   if (t < tce)
       phi = aphi * sin(2*pi*t/ tce) - bphi * sin(2*pi*t/ tce);
   else
       phi = 0;
   end

elastance = Emin * (1-phi) + Emax * phi;
end