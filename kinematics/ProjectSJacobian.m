function [ds_dqb,ds_dqc] = ProjectSJacobian(Linkage,q_b,q_c,s)

    T_r = FwdKinematicsAtS(Linkage,q_b,s);

    p_r = T_r(1:3,4);
    R_r = T_r(1:3,1:3);

    T_c = variable_expmap_g(q_c);
    p_c = T_c(1:3,4);

    r = p_r - p_c;

    t_full = dfzcds(Linkage,q_b,s);
    t = t_full(4:6);

    gs = dot(t,t);

    J_r = JacobianAtS(Linkage,q_b,s);

    Jp_r = R_r * J_r(4:6,:);

    g_qb = t.' * Jp_r;

    [~,Tg] = variable_expmap_gTg(q_c);

    Jp_c = T_c(1:3,1:3) * Tg(4:6,:);

    g_qc = -t.' * Jp_c;

    ds_dqb = -g_qb/gs;

    ds_dqc = -g_qc/gs;

end