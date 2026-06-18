function xi_s = xi_s_FD(Linkage, q_b,s)
    
    eps_fd = 1e-9;

    T_s = FwdKinematicsAtS(Linkage, q_b, s);
    T_p = FwdKinematicsAtS(Linkage, q_b, s+eps_fd);
    T_m = FwdKinematicsAtS(Linkage, q_b, s-eps_fd);

    xi_s = piecewise_logmap( ginv(T_s) * T_p ) - piecewise_logmap( ginv(T_s) * T_m );
    xi_s = xi_s / (2*eps);
end