function [t, qq1_MV, qq1_AO, qq1_PV] = interpolation_flux(whichQ)
if(whichQ == "piecewise_quadratic")
    t1 = linspace(0,0.2,200);
    t2 = linspace(0.2,0.45,250);
    t3 = linspace(0.45,0.50,50);
    t4 = linspace(0.50,0.55,50);
    t5 = linspace(0.55,0.60,50);
    t6 = linspace(0.60,1,400);
    tvec=[t1 t2 t3 t4 t5 t6];
    y1 = zeros(1,200);
    y2 =-18500*t2.*t2 + 12325*t2 - 1725;
    y3=-166666.666667*t3.*t3 + 156833.333333*t3-36750;
    y4=zeros(1,50);
    y5=-83333.3333333*t5.*t5+95833.3333333*t5-27500;
    y6=zeros(1,400);
    t=[t1 t2 t3 t4 t5 t6];
    qq1_MV=[y1 y2 y3 y4 y5 y6];
    dt=1/1000;
    qq1_AO = zeros(1000,1);
    qq1_PV = zeros(1000,1);

elseif(whichQ == "data_broome")
    load("mv_flow_broome")
    load("pv_flow_broome")
    load("ao_flow_broome")
    m3_to_ml=1000000;
    
    t_MV=mv_flow_broome(:,1);
    q_MV=mv_flow_broome(:,2);
    t_PV=pv_flow_broome(:,1);
    q_PV=pv_flow_broome(:,2);
    t_AO=ao_flow_broome(:,1);
    q_AO=ao_flow_broome(:,2);
    xq=linspace(0,1,1000);
    t=xq;
%     figure
    qq1_MV = interp1(t_MV,q_MV,xq);
%     subplot(3,1,1)
%     plot(t_MV,q_MV,'o',xq,qq1_MV,':.');
%     title("MV inflow")
%     xlabel("t(s)")
%     ylabel("Q (ml/s)")
%     
%     
%     % total_volume_MV=0;
%     % for(i=2:1000)
%     %     total_volume_MV= total_volume_MV + abs(qq_MV(i))*1/1000;
%     % end
%     % total_volume_MV
%     
%     subplot(3,1,2)
    qq1_PV = interp1(t_PV,q_PV,xq);
%     plot(t_PV,q_PV,'o',xq,qq1_PV,':.');
%     title("PV inflow")
%     xlabel("t(s)")
%     ylabel("Q (ml/s)")
%     total_volume_PV=0;
%     % for(i=2:1000)
%     %     total_volume_PV= total_volume_PV + abs(qq_PV(i))*1/1000;
%     % end
%     % total_volume_PV
%     
%     subplot(3,1,3)
    qq1_AO = interp1(t_AO,q_AO,xq);
%     plot(t_AO,q_AO,'o',xq,qq1_AO,':.');
%     title("AO inflow")
%     xlabel("t(s)")
%     ylabel("Q (ml/s)")
elseif(whichQ == "data_broome_MVopen_active")
    load("mv_flow_broome")
    load("pv_flow_broome")
    load("ao_flow_broome")
    m3_to_ml=1000000;
    
    t_MV=mv_flow_broome(:,1);
    q_MV=mv_flow_broome(:,2);
    t_PV=pv_flow_broome(:,1);
    q_PV=pv_flow_broome(:,2);
    t_AO=ao_flow_broome(:,1);
    q_AO=ao_flow_broome(:,2);
    xq=linspace(0,1,1000);
    t=xq;
    qq1_MV_aux = interp1(t_MV,q_MV,xq);
    qq1_PV_aux = interp1(t_PV,q_PV,xq);
    qq1_AO_aux = interp1(t_AO,q_AO,xq);
    qq1_MV = zeros(1,1000);
    qq1_PV = zeros(1,1000);
    qq1_AO = zeros(1,1000);
    qq1_MV(1,1:351) = qq1_MV_aux(650:1000);
    qq1_MV(1,352:1000) = qq1_MV_aux(1:649);
    qq1_PV(1,1:351) = qq1_PV_aux(650:1000);
    qq1_PV(1,352:1000) = qq1_PV_aux(1:649);
    qq1_AO(1,1:351) = qq1_AO_aux(650:1000);
    qq1_AO(1,352:1000) = qq1_AO_aux(1:649);
    elseif(whichQ == "data_broome_MVopen_passive")
    load("mv_flow_broome")
    load("pv_flow_broome")
    load("ao_flow_broome")
    m3_to_ml=1000000;
    
    t_MV=mv_flow_broome(:,1);
    q_MV=mv_flow_broome(:,2);
    t_PV=pv_flow_broome(:,1);
    q_PV=pv_flow_broome(:,2);
    t_AO=ao_flow_broome(:,1);
    q_AO=ao_flow_broome(:,2);
    xq=linspace(0,1,1000);
    t=xq;
    qq1_MV_aux = interp1(t_MV,q_MV,xq);
    qq1_PV_aux = interp1(t_PV,q_PV,xq);
    qq1_AO_aux = interp1(t_AO,q_AO,xq);
    qq1_MV = zeros(1,1000);
    qq1_PV = zeros(1,1000);
    qq1_AO = zeros(1,1000);
    qq1_MV(1,1:731) = qq1_MV_aux(1,270:1000);
    qq1_MV(1,732:1000) = qq1_MV_aux(1,1:269);
    qq1_PV(1,1:731) = qq1_PV_aux(1,270:1000);
    qq1_PV(1,732:1000) = qq1_PV_aux(1,1:269);
    qq1_AO(1,1:731) = qq1_AO_aux(1,270:1000);
    qq1_AO(1,732:1000) = qq1_AO_aux(1,1:269);
elseif(whichQ == "data_prisco")
    load("mv_flow_prisco")   
    t_MV=mv_flow_prisco(:,1);
    q_MV=mv_flow_prisco(:,2);
    xq=linspace(0,1,1000);
    t=xq;
    qq1_MV = interp1(t_MV,q_MV,xq);
    qq1_AO = zeros(1000,1);
    qq1_PV = zeros(1000,1);
    
else
    fprintf("No valid option for prescribed flux");
end
end