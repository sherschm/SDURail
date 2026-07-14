function dqse3dqm = Screw2EulerJacobian(q)

    %find q_se3
    %theta = q(4:6);
    r = q(1:3);
    Rm = computeRotationMatrix(q);
    T_floating = [Rm r;0 0 0 1];
    
    q_se3 = piecewise_logmap(T_floating);

    Jl_inv = inv(SO3LeftJacobian(q_se3(1:3)));
    drho_se3_dr = Jl_inv;
    
    dtheta_se3_dtheta = [LogSO3(Rm'*computeRotationMatrixDerivative(q, 4)), ...
        LogSO3(Rm'*computeRotationMatrixDerivative(q, 5)), ...
        LogSO3(Rm'*computeRotationMatrixDerivative(q, 6))];

    dqse3dqm = [zeros(3,3) dtheta_se3_dtheta; ...
                drho_se3_dr zeros(3,3)];
end