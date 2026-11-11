% 
% function [e_out,J_out,Jd_out] = ErrorJAtS(Linkage, Carriage, s, q, qd) 
%     S = [1 0 0 0 0 0; ...
%         0 1 0 0 0 0; ...
%         0 0 1 0 0 0; ...
%         0 0 0 0 1 0; ...
%         0 0 0 0 0 1];
% 
%     %S= eye(6);
% 
%     n_beam = Linkage.ndof;
%     n_mass = Carriage.ndof;
%     q_b = q(1:n_beam);
%     q_mass = q(n_beam+1:n_beam+n_mass);
%     qd_b = qd(1:n_beam);
%     qd_mass = qd(n_beam+1:n_beam+n_mass);
% 
%     % --- Forward kinematics 
%     T_s = FwdKinematicsAtS(Linkage, q_b, s);
%     [T_m,Jm_WRONG,Jmd_WRONG] = variable_expmap_gTgTgd(q_mass,qd_mass);
% 
%     Jm = zeros(6,6);
%     Jm(1:3,1:3) = Jm_WRONG(1:3,1:3)';
%     Jm(4:6,4:6) = Jm_WRONG(4:6,4:6)';
%     Jm(4:6,1:3) = Jm_WRONG(4:6,1:3)';
% 
%     Jmd = zeros(6,6);
%     Jmd(1:3,1:3) = Jmd_WRONG(1:3,1:3)';
%     Jmd(4:6,4:6) = Jmd_WRONG(4:6,4:6)';
%     Jmd(4:6,1:3) = Jmd_WRONG(4:6,1:3)';
% 
%     %T_m = variable_expmap_g(q_mass); 
%     g_sm = ginv(T_s) * T_m;
%     %g_sm = ginv(T_m) * T_s;
% 
%     % --- Jacobians in their natural frames
%     Js = JacobianAtS(Linkage, q_b, s);
%     Jsd = JacobiandotAtS(Linkage, q_b,qd_b, s);
%     %J_m = SE3RightJacobianFromPose(q_mass); 
% 
%     % --- Adjoint terms
%     Ad_sm_inv = dinamico_Adjoint(ginv(g_sm)); 
% 
%     %T_s = FwdKinematicsAtS(Linkage,q_b,s);
%     %T_m = variable_expmap_g(q_mass);
% 
%     e  = piecewise_logmap(g_sm);
%     Jr = SE3RightJacobianFromPose(e);
%     Jr_inv = Jr \ eye(6);
% 
%     %Jb =  J_s ;
%     %Jm =  - Ad_sm_inv *J_m;
%     %Jm =  - Ad_sm_inv *J_m;
%     %Jm =  - J_m;
% 
%     %dphi_ds = dphi_ds_FD(Linkage,q_b,q_mass,s);
%     %xi_s = xi_s_FD(Linkage, q_b,s);
%     %[xi_s,dxi_ds] = dfzcds(Linkage,q,s);
% 
%     %dphi_ds =  invLeftJacobianSE3(e)*Ad_sm_inv*xi_s; %dphi_ds_FD_Lie(Linkage, q_b, q_mass, s);
%     %[ds_dqb, ds_dqm] = ProjectSJacobianFD(Linkage,q_b,q_mass,s);
%     %J = Jr_inv * [-Ad_sm_inv*(J_s+dphi_ds*ds_dqb) (Jm-Ad_sm_inv*dphi_ds*ds_dqm)];
% 
%     %J =  Jr_inv*[-Ad_sm_inv*Js Jm];
%     J =  [-Js Jm];
% 
%     %J =  [-Js Jm];
% 
%     Jd =  [-Jsd Jmd];
%     %Jd =  zeros(size(J));
% 
%     e_out = S*e;
%     J_out = S*J;
%     Jd_out = S*Jd;
% 
% 
% end


function [e_out,J_out,Jd_out] = ErrorJAtS(Linkage, Carriage, s, q, qd, opts) 
arguments
    Linkage
    Carriage
    s
    q
    qd
    opts.fixed logical = false
end

if opts.fixed
    S = eye(6);
else
    S = [1 0 0 0 0 0; ...
        0 1 0 0 0 0; ...
        0 0 1 0 0 0; ...
        0 0 0 0 1 0; ...
        0 0 0 0 0 1];
end

%S= eye(6);

T_offset = Carriage.VLinks(1).gi;

n_beam = Linkage.ndof;
n_mass = Carriage.ndof;
q_b = q(1:n_beam);
q_mass = q(n_beam+1:n_beam+n_mass);
qd_b = qd(1:n_beam);
qd_mass = qd(n_beam+1:n_beam+n_mass);

% --- Forward kinematics 
T_s = FwdKinematicsAtS(Linkage, q_b, s);

[T_m,Jm_WRONG,Jmd_WRONG] = variable_expmap_gTgTgd(q_mass,qd_mass);

T_m_offset = T_m*T_offset;
Ad_offset = dinamico_Adjoint(ginv(T_offset));

Jm_WRONG = Ad_offset*Jm_WRONG;
Jmd_WRONG = Ad_offset*Jmd_WRONG;

Jm = zeros(6,6);
Jm(1:3,1:3) = Jm_WRONG(1:3,1:3)';
Jm(4:6,4:6) = Jm_WRONG(4:6,4:6)';
Jm(4:6,1:3) = Jm_WRONG(4:6,1:3)';

Jmd = zeros(6,6);
Jmd(1:3,1:3) = Jmd_WRONG(1:3,1:3)';
Jmd(4:6,4:6) = Jmd_WRONG(4:6,4:6)';
Jmd(4:6,1:3) = Jmd_WRONG(4:6,1:3)';

g_sm = ginv(T_s) * T_m_offset;

% --- Jacobians in their natural frames
Js = JacobianAtS(Linkage, q_b, s);
Jsd = JacobiandotAtS(Linkage, q_b,qd_b, s);

e  = piecewise_logmap(g_sm);
J =  [-Js Jm];
Jd =  [-Jsd Jmd];

e_out = S*e;
J_out = S*J;
Jd_out = S*Jd;

end
