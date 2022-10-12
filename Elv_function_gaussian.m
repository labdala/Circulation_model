function elastance = Elv_function_gaussian(t, a , b, c, E_lv_base)
    elastance = 1/(a * exp(-(t-b)^2/(2*c^2)) + E_lv_base);
end