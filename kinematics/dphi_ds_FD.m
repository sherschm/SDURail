function dphi_ds = dphi_ds_FD(Linkage, q_b,q_mass,s)
    
    nc = 5;

    eps_fd = 1e-4;
    T_s = FwdKinematicsAtS(Linkage, q_b, s+eps_fd);
    T_m = variable_expmap_g(q_mass); 
    g_sm = ginv(T_s) * T_m;
    phi_p = piecewise_logmap(g_sm);

    T_s = FwdKinematicsAtS(Linkage, q_b, s-eps_fd);
    T_m = variable_expmap_g(q_mass); 
    g_sm = ginv(T_s) * T_m;
    phi_m = piecewise_logmap(g_sm);

    dphi_ds = (phi_p-phi_m)/(2*eps_fd);

end