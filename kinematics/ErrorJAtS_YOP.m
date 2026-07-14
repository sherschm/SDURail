
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
    
    %convert Euler coordinates to se(3)


    % --- Forward kinematics 
    T_s = FwdKinematicsAtS(Linkage, q_b, s);
    T_m = [computeRotationMatrix(q_mover) q_mover(1:3);0 0 0 1];
    
    q_mover_se3 = piecewise_logmap(T_m);
    qd_mover_se3 = zeros(6); %CORRECT THIS!!!

    [foo,Jm,Jmd] = variable_expmap_gTgTgd(q_mover_se3,qd_mover_se3);
    
    %T_m = variable_expmap_g(q_mover); 
    g_sm = ginv(T_s) * T_m;

    % --- Jacobians in their natural frames
    Js = JacobianAtS(Linkage, q_b, s);
    Jsd = JacobiandotAtS(Linkage, q_b,qd_b, s);

    % --- Adjoint terms
    Ad_sm_inv = dinamico_Adjoint(ginv(g_sm)); 

    e  = piecewise_logmap(g_sm);
    Jr = SE3RightJacobianFromPose(e);
    Jr_inv = Jr \ eye(6);

    J =  [-Js Jm*J_se3_to_euler];

    Jd =  [-Jsd Jmd*zeros(size(J_se3_to_euler))]; %CORRECT THIS!

    e_out = S*e;
    J_out = S*J;
    Jd_out = S*Jd;
end
