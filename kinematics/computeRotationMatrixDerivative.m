function dA = computeRotationMatrixDerivative(Z, idx)

phi   = Z(4);
theta = Z(5);
psi   = Z(6);

Rx = [1, 0, 0;
      0, cos(phi), -sin(phi);
      0, sin(phi),  cos(phi)];

Ry = [ cos(theta), 0, sin(theta);
       0,          1, 0;
      -sin(theta), 0, cos(theta)];

Rz = [cos(psi), -sin(psi), 0;
      sin(psi),  cos(psi), 0;
      0,         0,        1];

dRx = [0, 0, 0;
       0, -sin(phi), -cos(phi);
       0,  cos(phi), -sin(phi)];

dRy = [-sin(theta), 0, cos(theta);
        0,          0, 0;
       -cos(theta), 0, -sin(theta)];

dRz = [-sin(psi), -cos(psi), 0;
        cos(psi), -sin(psi), 0;
        0,         0,        0];

switch idx
    case 4      % dA/dphi
        dA = Rz * Ry * dRx;

    case 5      % dA/dtheta
        dA = Rz * dRy * Rx;

    case 6      % dA/dpsi
        dA = dRz * Ry * Rx;

    otherwise
        error('idx must be 4, 5, or 6.');
end

end