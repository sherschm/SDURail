function Jrd = SE3RightJacobianDot_FD(xi, xidot)

eps = 1e-8;

J0 = SE3RightJacobianFromPose(xi);
J1 = SE3RightJacobianFromPose(xi + eps*xidot);

Jrd = (J1 - J0)/eps;

end