
function [e_out,J_out,Jd_out] = ErrorJAtS_YOP(Linkage, Carriage, s, q, qd) 
    S = [1 0 0 0 0 0; ...
        0 1 0 0 0 0; ...
        0 0 1 0 0 0; ...
        0 0 0 0 1 0; ...
        0 0 0 0 0 1];

    %S= eye(6);

    n_beam = Linkage.ndof;
    n_mover = 7; %CORRECT THIS
    q_b = q(1:n_beam);
    q_mover = q(n_beam+1:n_beam+n_mover);

    qd_b = qd(1:n_beam);
    qd_mover = qd(n_beam+1:n_beam+n_mover);
    

    % --- Forward kinematics 
    T_s = FwdKinematicsAtS(Linkage, q_b, s);
    r_m = q_mover(1:3);
    T_m = [computeRotationMatrix(q_mover) r_m(:);0 0 0 1];
    
    % --- Convert to se(3)
    q_mover_se3 = piecewise_logmap(T_m);
    qd_mover_se3 = zeros(6); %CORRECT THIS!!!

    [foo,Jm,Jmd] = variable_expmap_gTgTgd(q_mover_se3,qd_mover_se3);
    
    g_sm = ginv(T_s) * T_m;

    % --- Jacobians 
    Js = JacobianAtS(Linkage, q_b, s);
    Jsd = JacobiandotAtS(Linkage, q_b,qd_b, s);

    % --- Adjoint terms
    Ad_sm_inv = dinamico_Adjoint(ginv(g_sm)); 
    
    e  = piecewise_logmap(g_sm);
    Jr = SE3RightJacobianFromPose(e);
    Jr_inv = Jr \ eye(6);
    
    J_se3_to_euler = [Screw2EulerJacobian(q_mover) zeros(6,n_mover-6)];

    %J =  Jr_inv*[-Ad_sm_inv*Js Jm];
    J =  [-Js Jm*J_se3_to_euler];

    %J =  [-Js Jm];

    Jd =  [-Jsd Jmd*zeros(size(J_se3_to_euler))]; %CORRECT THIS!
    %Jd =  zeros(size(J));

    e_out = S*e;
    J_out = S*J;
    Jd_out = S*Jd;
end
