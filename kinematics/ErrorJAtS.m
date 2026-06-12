
function [e_out,J_out] = ErrorJAtS(Linkage, Carriage, s, q, xi_rs) 

    n_beam = Linkage.ndof;
    n_mass = 6;
    q_b = q(1:n_beam);
    q_mass = q(n_beam+1:n_beam+n_mass);
    %qd_b = qd(1:n_beam);
    %qd_mass = qd(n_beam+1:n_beam+n_mass);
    
    S = [eye(3) zeros(3,3);zeros(2,4) eye(2)];

    % --- Forward kinematics 
    T_s = FwdKinematicsAtS(Linkage, q_b, s);
    T_m = variable_expmap_g(q_mass);
    g_sm = ginv(T_s) * T_m;
    
    T_rs_inv = ginv(variable_expmap_g(xi_rs));
    dT1 = T_rs_inv*T_s;
    dT2 = T_rs_inv*T_m;
    
    dxi1 = piecewise_logmap(dT1);
    dxi2 = piecewise_logmap(dT2);
    
    e_out = [dxi1;S*dxi2];
    
    Jr_xi1_inv = SE3RightJacobianFromPose(dxi1)\eye(6);
    Jr_xi2_inv = SE3RightJacobianFromPose(dxi2)\eye(6);
    
    dphim1_dqr = Jr_xi1_inv * JacobianAtS(Linkage, q_b, s);
    dphim1_ds = Jr_xi1_inv * dfzcds(Linkage,q_b,s);
    dphim1_dxirs = - Jr_xi1_inv;

    dphim2_dqm = S* Jr_xi2_inv;
    dphim2_dxirs = -S * Jr_xi2_inv*dinamico_Adjoint(ginv(dT2));
    

    dphim1_dx = [dphim1_dqr zeros(6,6) dphim1_ds dphim1_dxirs];
    dphim2_dx = [zeros(5,n_beam) dphim2_dqm zeros(5,1) dphim2_dxirs];

    J_out = [dphim1_dx;dphim2_dx];
   
end

% function [e_out,J_out,Jd_out] = ErrorJAtS(Linkage, Carriage, s, q, qd)
% 
%     n_beam = Linkage.ndof;
%     n_mass = 6;
% 
%     q_b    = q(1:n_beam);
%     q_mass = q(n_beam+1:n_beam+n_mass);
% 
%     qd_b    = qd(1:n_beam);
%     qd_mass = qd(n_beam+1:n_beam+n_mass);
% 
%     % --- Forward kinematics
%     T_s = FwdKinematicsAtS(Linkage, q_b, s);
%     T_m = variable_expmap_g(q_mass);
% 
%     g_sm = ginv(T_s) * T_m;
% 
%     % =========================================================
%     % ERROR (same as before)
%     % =========================================================
%     e = piecewise_logmap(g_sm);
% 
%     % =========================================================
%     % BODY TWISTS (consistent with closed-loop formulation)
%     % =========================================================
%     J_s  = JacobianAtS(Linkage, q_b, s);
%     Jd_s = JacobiandotAtS(Linkage, q_b, qd_b, s);
% 
%     J_m  = SE3RightJacobianFromPose(q_mass);
%     %disp(J_m)
%     Jd_m = SE3RightJacobianDot(q_mass, qd_mass);
% 
%     % --- body velocities
%     Vs = J_s * qd_b;
%     Vm = J_m * qd_mass;
% 
%     % =========================================================
%     % RELATIVE VELOCITY (THIS IS THE KEY CONSISTENT FORM)
%     % =========================================================
%     Ad_sm_inv = dinamico_Adjoint(ginv(g_sm));
%     Ad_s_inv  = dinamico_Adjoint(ginv(T_s));
% 
%     Vrel = -Ad_sm_inv * Vs + Ad_s_inv * Vm;
% 
%     % =========================================================
%     % JACOBIAN (NO J_log)
%     % =========================================================
%     J = [-Ad_sm_inv * J_s , Ad_s_inv * J_m];
% 
%     % =========================================================
%     % TIME DERIVATIVE OF JACOBIAN (CONSISTENT VERSION)
%     % =========================================================
% 
%     Vs_dot = Jd_s * qd_b;   % (approx consistent with your code)
%     Vm_dot = Jd_m * qd_mass;
% 
%     % adjoint time derivatives (same structure as closed-loop code logic)
%     Ad_sm_inv_dot = -Ad_sm_inv * dinamico_adj(Vrel);
%     Ad_s_inv_dot  = -Ad_s_inv  * dinamico_adj(Vs);
% 
%     Jrel_dot = [ ...
%         -(Ad_sm_inv_dot * J_s + Ad_sm_inv * Jd_s), ...
%           (Ad_s_inv_dot  * J_m + Ad_s_inv  * Jd_m)
%     ];
% 
%     %Jd = Jrel_dot;
%     Jd = zeros(6,n_beam+n_mass);
% 
%    J_out = [J(1:3,:) ;J(5:6,:)];
%    Jd_out = [Jd(1:3,:) ;Jd(5:6,:)];
%    e_out = [e(1:3); e(5:6)];
% 
% end