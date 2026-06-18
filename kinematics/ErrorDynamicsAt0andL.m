function [err0,J0,Jd0,errL,JL,JdL] = ErrorDynamicsAt0andL( ...
    Linkage, Carriage, q_b, qd_b, T_0_fixed, T_L_fixed)

    % ----------------------------
    % Kinematics
    % ----------------------------
    T_full = FwdKinematics(Linkage, q_b);

    T_0 = T_full(1:4,:);
    T_L = T_full(end-3:end,:);

    % ----------------------------
    % Consistent BODY-frame errors
    % e = log( T^{-1} * T_des )
    % ----------------------------
    %dT_0 = ginv(T_0) * T_0_fixed;
    %dT_L = ginv(T_L) * T_L_fixed;
    dT_0 = ginv(T_0_fixed) * T_0;
    dT_L = ginv(T_L_fixed) * T_L;

    err0 = piecewise_logmap(dT_0);
    errL = piecewise_logmap(dT_L);

    % ----------------------------
    % Body Jacobians (already returned)
    % ----------------------------
    J_full = Jacobian(Linkage, q_b);

    J0 = [J_full(1:6,:), zeros(6, Carriage.ndof)];
    JL = [J_full(end-5:end,:), zeros(6, Carriage.ndof)];

    % ----------------------------
    % Time derivative of BODY Jacobian
    % Must be consistent with Jacobian()
    % ----------------------------
    Jd_full = Jacobiandot_corrected(Linkage, q_b, qd_b);

    Jd0 = [Jd_full(1:6,:), zeros(6, Carriage.ndof)];
    JdL = [Jd_full(end-5:end,:), zeros(6, Carriage.ndof)];

    % spherical Pin-spherical Pin
    % S = [zeros(3,3) eye(3)];
    % err0 = S*err0;
    % J0 = S*J0;
    % Jd0 = S*Jd0;
    % errL = S*errL;
    % JL = S*JL;
    % JdL = S*JdL;

    % Fixed-Fixed
    S = eye(6);
    err0 = S*err0;
    J0 = S*J0;
    Jd0 = S*Jd0;
    errL = S*errL;
    JL = S*JL;
    JdL = S*JdL;
end