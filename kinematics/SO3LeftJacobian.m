function J = SO3LeftJacobian(phi)
%LEFTJACOBIANSO3 Computes the SO(3) left Jacobian.
%
% Input:
%   phi : 3x1 rotation vector
%
% Output:
%   J   : 3x3 left Jacobian of SO(3)

theta = norm(phi);

Phi = dinamico_tilde(phi);

% Small-angle approximation
if theta < 1e-8
    % Taylor expansion
    J = eye(3) ...
        + 0.5 * Phi ...
        + (1/6) * (Phi * Phi);
else
    A = (1 - cos(theta)) / theta^2;
    B = (theta - sin(theta)) / theta^3;

    J = eye(3) ...
        + A * Phi ...
        + B * (Phi * Phi);
end

end
