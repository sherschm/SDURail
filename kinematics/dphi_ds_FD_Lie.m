function dphi_ds = dphi_ds_FD_Lie(Linkage, q_b, q_mass, s)

    eps_fd = 1e-6;  % smaller is usually better for Lie FD

    % --- central points ---
    Tm = variable_expmap_g(q_mass);

    % =========================
    % + perturbation
    % =========================
    T_s_p = FwdKinematicsAtS(Linkage, q_b, s + eps_fd);
    g_p = ginv(T_s_p) * Tm;
    phi_p = piecewise_logmap(g_p);

    % =========================
    % - perturbation
    % =========================
    T_s_m = FwdKinematicsAtS(Linkage, q_b, s - eps_fd);
    g_m = ginv(T_s_m) * Tm;
    phi_m = piecewise_logmap(g_m);

    % =========================
    % central difference in Lie algebra
    % =========================
    dphi_raw = (phi_p - phi_m) / (2 * eps_fd);

    % =========================
    % Lie correction (important)
    % =========================
    phi0 = piecewise_logmap(ginv(FwdKinematicsAtS(Linkage, q_b, s)) * Tm);

    Jl_inv = invLeftJacobianSE3(phi0);

    dphi_ds = Jl_inv * dphi_raw;

end
