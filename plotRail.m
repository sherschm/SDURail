
function plotRail(Linkage,x)
    
    plotq(Linkage,x(1:Linkage.ndof));
    hold on;
    g_s = variable_expmap_g(x(Linkage.ndof+1:Linkage.ndof+6));
    %FwdKinematicsAtS(RailLinkage,q_b0,s0);
    plotFrame(g_s)
    drawnow;
end
