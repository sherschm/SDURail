
function [e_out,J_out] = ErrorJAtS(Linkage, Carriage, s, q) 
    S = [1 0 0 0 0 0; ...
        0 1 0 0 0 0; ...
        0 0 1 0 0 0; ...
        0 0 0 0 1 0; ...
        0 0 0 0 0 1];

    n_beam = Linkage.ndof;
    n_mass = Carriage.ndof;
    q_b = q(1:n_beam);
    q_mass = q(n_beam+1:n_beam+n_mass);

    % --- Forward kinematics 
    T_s = FwdKinematicsAtS(Linkage, q_b, s);
    T_m = variable_expmap_g(q_mass); 
    g_sm = ginv(T_s) * T_m;

    % --- Jacobians in their natural frames
    J_s = JacobianAtS(Linkage, q_b, s);
    J_m = SE3RightJacobianFromPose(q_mass); 

    % --- Adjoint terms
    Ad_sm_inv = dinamico_Adjoint(ginv(g_sm)); 

    T_s = FwdKinematicsAtS(Linkage,q_b,s);
    T_m = variable_expmap_g(q_mass);
    
    e  = piecewise_logmap(ginv(T_s)*T_m);
    Jr = SE3RightJacobianFromPose(e);
    Jr_inv = Jr \ eye(6);
    
    Jb = - Ad_sm_inv * J_s ;
    Jm =  J_m;

    %dphi_ds = dphi_ds_FD(Linkage,q_b,q_mass,s);
    %xi_s = xi_s_FD(Linkage, q_b,s);
    %[xi_s,dxi_ds] = dfzcds(Linkage,q,s);

    %dphi_ds =  invLeftJacobianSE3(e)*Ad_sm_inv*xi_s; %dphi_ds_FD_Lie(Linkage, q_b, q_mass, s);
    %[ds_dqb, ds_dqm] = ProjectSJacobianFD(Linkage,q_b,q_mass,s);
    %J = Jr_inv * [-Ad_sm_inv*(J_s+dphi_ds*ds_dqb) (Jm-Ad_sm_inv*dphi_ds*ds_dqm)];
    J = Jr_inv * [Jb Jm];

    J_out = S*J;
    e_out = S*e;
end
