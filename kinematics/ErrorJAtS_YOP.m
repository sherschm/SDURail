
function [e_out,J_out,Jd_out] = ErrorJAtS(Linkage, Carriage, s, q, qd) 
    S = [1 0 0 0 0 0; ...
        0 1 0 0 0 0; ...
        0 0 1 0 0 0; ...
        0 0 0 0 1 0; ...
        0 0 0 0 0 1];

    %S= eye(6);

    n_beam = Linkage.ndof;
    n_mass = Carriage.ndof;
    q_b = q(1:n_beam);
    q_mover = q(n_beam+1:n_beam+n_mass);
    qd_b = qd(1:n_beam);
    qd_mass = qd(n_beam+1:n_beam+n_mass);

    % --- Forward kinematics 
    T_s = FwdKinematicsAtS(Linkage, q_b, s);
    [T_m,Jm,Jmd] = variable_expmap_gTgTgd(q_mover,qd_mass);

    %T_m = variable_expmap_g(q_mover); 
    g_sm = ginv(T_s) * T_m;
    %g_sm = ginv(T_m) * T_s;

    % --- Jacobians in their natural frames
    Js = JacobianAtS(Linkage, q_b, s);
    Jsd = JacobiandotAtS(Linkage, q_b,qd_b, s);
    %J_m = SE3RightJacobianFromPose(q_mover); 

    % --- Adjoint terms
    Ad_sm_inv = dinamico_Adjoint(ginv(g_sm)); 

    %T_s = FwdKinematicsAtS(Linkage,q_b,s);
    %T_m = variable_expmap_g(q_mover);
    
    e  = piecewise_logmap(g_sm);
    Jr = SE3RightJacobianFromPose(e);
    Jr_inv = Jr \ eye(6);
    
    %Jb =  J_s ;
    %Jm =  - Ad_sm_inv *J_m;
    %Jm =  - Ad_sm_inv *J_m;
    %Jm =  - J_m;

    %dphi_ds = dphi_ds_FD(Linkage,q_b,q_mover,s);
    %xi_s = xi_s_FD(Linkage, q_b,s);
    %[xi_s,dxi_ds] = dfzcds(Linkage,q,s);

    %dphi_ds =  invLeftJacobianSE3(e)*Ad_sm_inv*xi_s; %dphi_ds_FD_Lie(Linkage, q_b, q_mover, s);
    %[ds_dqb, ds_dqm] = ProjectSJacobianFD(Linkage,q_b,q_mover,s);
    %J = Jr_inv * [-Ad_sm_inv*(J_s+dphi_ds*ds_dqb) (Jm-Ad_sm_inv*dphi_ds*ds_dqm)];
    
    J_se3_to_euler = [Screw2EulerJacobian(q_mover) zeros(6,1)];

    %J =  Jr_inv*[-Ad_sm_inv*Js Jm];
    J =  [-Js Jm*J_se3_to_euler];

    %J =  [-Js Jm];

    Jd =  [-Jsd Jmd*zeros(size(J_se3_to_euler))];
    %Jd =  zeros(size(J));

    e_out = S*e;
    J_out = S*J;
    Jd_out = S*Jd;
end
