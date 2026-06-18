
function Jinv = invLeftJacobianSE3(xi)

    % xi = [rho; phi] in se(3)

    rho = xi(1:3);
    phi = xi(4:6);

    theta = norm(phi);

    if theta < 1e-8
        Jinv = eye(6);
        return;
    end

    % skew operator
    W = dinamico_tilde(phi);

    A = (theta/2) * cot(theta/2);
    B = 1 - A;

    Jinv_so3 = eye(3) + 0.5*W + (1/theta^2)*(1 - A/2)*W*W;

    V = dinamico_tilde(rho);

    Jinv = [Jinv_so3, zeros(3);
            0.5*V + (1/theta^2)*(1 - A/2)*(W*V + V*W), Jinv_so3];

end