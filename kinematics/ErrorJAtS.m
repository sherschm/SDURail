function [e,J,Jd] = ErrorJAtS(Linkage, Carriage, s, q, qd)

n_beam = Linkage.ndof;
n_mass = 6;

q_b    = q(1:n_beam);
q_mass = q(n_beam+1:n_beam+n_mass);

qd_b    = qd(1:n_beam);
qd_mass = qd(n_beam+1:n_beam+n_mass);

% --- Forward kinematics
T_s = FwdKinematicsAtS(Linkage, q_b, s);
T_m = variable_expmap_g(q_mass);

g_sm = ginv(T_s) * T_m;

% =========================================================
% ERROR (same as before)
% =========================================================
e = piecewise_logmap(g_sm);

% =========================================================
% BODY TWISTS (consistent with closed-loop formulation)
% =========================================================
J_s  = JacobianAtS(Linkage, q_b, s);
Jd_s = JacobiandotAtS(Linkage, q_b, qd_b, s);

J_m  = SE3RightJacobianFromPose(q_mass);
Jd_m = SE3RightJacobianDot(q_mass, qd_mass);

% --- body velocities
Vs = J_s * qd_b;
Vm = J_m * qd_mass;

% =========================================================
% RELATIVE VELOCITY (THIS IS THE KEY CONSISTENT FORM)
% =========================================================
Ad_sm_inv = dinamico_Adjoint(ginv(g_sm));
Ad_s_inv  = dinamico_Adjoint(ginv(T_s));

Vrel = -Ad_sm_inv * Vs + Ad_s_inv * Vm;

% =========================================================
% JACOBIAN (NO J_log)
% =========================================================
J = [-Ad_sm_inv * J_s , Ad_s_inv * J_m];

% =========================================================
% TIME DERIVATIVE OF JACOBIAN (CONSISTENT VERSION)
% =========================================================

Vs_dot = Jd_s * qd_b;   % (approx consistent with your code)
Vm_dot = Jd_m * qd_mass;

% adjoint time derivatives (same structure as closed-loop code logic)
Ad_sm_inv_dot = -Ad_sm_inv * dinamico_adj(Vrel);
Ad_s_inv_dot  = -Ad_s_inv  * dinamico_adj(Vs);

Jrel_dot = [ ...
    -(Ad_sm_inv_dot * J_s + Ad_sm_inv * Jd_s), ...
      (Ad_s_inv_dot  * J_m + Ad_s_inv  * Jd_m)
];

Jd = Jrel_dot;

end