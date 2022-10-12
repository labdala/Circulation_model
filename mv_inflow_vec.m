function [t, d_Q_mv] = mv_inflow_vec(which_Q_mv, ncycles, n_points, period, t0)
    if (which_Q_mv == "mv_flow_prisco")
        load("mv_flow_prisco")   
        t_MV1=(mv_flow_prisco(54:86,1)-mv_flow_prisco(54,1))*0.1/(mv_flow_prisco(86,1)-mv_flow_prisco(54,1));
        q_MV1=mv_flow_prisco(54:86,2);
        t_MV3_aux=[((mv_flow_prisco(94:98,1)-0.4783))' (mv_flow_prisco(2:53,1)+mv_flow_prisco(98,1)-0.4783)'];
        t_MV3 = t_MV3_aux(1,1) + (t_MV3_aux-t_MV3_aux(1,1)).*0.3./(t_MV3_aux(1,57)-t_MV3_aux(1,1));
        q_MV3=[(mv_flow_prisco(94:98,2))' (mv_flow_prisco(2:53,2))'];
        xq1=linspace(0,0.1,125);
        xq2=linspace(0.1,0.5,500);
        xq3=linspace(0.5,period,375);
        qq1_MV1 = interp1(t_MV1,q_MV1,xq1);
        qq1_MV2= 0*xq2;
        qq1_MV3 = interp1(t_MV3,q_MV3,xq3);
        qq1_MV3(375)=qq1_MV1(1);
        tt=[xq1 xq2 xq3];
        d_Q_mv0 = 1.3*[qq1_MV1 qq1_MV2 qq1_MV3];
    elseif(which_Q_mv == "constant")
        int = fix(n_points/ncycles);
        d_Q_mv0=100*ones(1,1000);
        tt = linspace(t0, period,1000 );
    else
        error("The selected input MV flow has not been implemented");
    end

t=tt;
d_Q_mv = d_Q_mv0;
for i = 1 : ncycles-1
    t = cat(2, t, tt + period * i);
    d_Q_mv = cat(2, d_Q_mv, d_Q_mv0); 
end
