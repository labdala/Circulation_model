function [t, pp1_MV,pp1_LV,pp1_AO] = interpolation_pressure(whichP)
if(whichP == "mmhg_pironet")
    load("la_pressure_mmhg_pironet")
    
    t_MV=la_pressure_mmhg_pironet(:,1);
    p_MV=la_pressure_mmhg_pironet(:,2);
    xq=linspace(0,1,1000);
    t=xq;
    pp1_MV = interp1(t_MV,p_MV,xq);
elseif(whichP=="mmhg_pironet_MVopen_passive")
    load("la_pressure_mmhg_pironet")
    
    t_MV=la_pressure_mmhg_pironet(:,1);
    p_MV=la_pressure_mmhg_pironet(:,2);
    xq=linspace(0,1,1000);
    t=xq;
    pp1_MV_aux = interp1(t_MV,p_MV,xq);
    pp1_MV(1,1:600) = pp1_MV_aux(401:1000);
    pp1_MV(1,601:1000) = pp1_MV_aux(1:400);
    pp1_LV=zeros(1000,1);
    pp1_AO=zeros(1000,1);
elseif(whichP=="mmhg_braunwnvald")
    load("braunwnvald_la_pressure")
    load("braunwnvald_ao_pressure")
    load("braunwnvald_lv_pressure")
    
    t_LA=braunwnvald_la_pressure(:,1);
    p_LA=braunwnvald_la_pressure(:,2);
    t_LV=braunwnvald_lv_pressure(:,1);
    p_LV=braunwnvald_lv_pressure(:,2);
    t_AO=braunwnvald_ao_pressure(:,1);
    p_AO=braunwnvald_ao_pressure(:,2);
    xq=linspace(0,3,1000);
    t=xq;
    pp1_LA = interp1(t_LA,p_LA,xq);
    pp1_LV = interp1(t_LV,p_LV,xq);
    pp1_AO = interp1(t_AO,p_AO,xq);
    pp1_MV = pp1_LA;
elseif(whichP=="mmhg_wiggers")    
    load("wiggers_la_pressure.csv")
    load("wiggers_ao_pressure.csv")
    load("wiggers_lv_pressure.csv")
    
    t_LA=wiggers_la_pressure(:,1);
    p_LA=wiggers_la_pressure(:,2);
    t_LV=wiggers_lv_pressure(:,1);
    p_LV=wiggers_lv_pressure(:,2);
    t_AO=wiggers_ao_pressure(:,1);
    p_AO=wiggers_ao_pressure(:,2);
    xq=linspace(0,0.8,1000);
    t=xq;
    pp1_LA = interp1(t_LA,p_LA,xq);
    pp1_LV = interp1(t_LV,p_LV,xq);
    pp1_AO = interp1(t_AO,p_AO,xq);
    pp1_MV = pp1_LA;
end
end