function dx = ODESolverImplicit(t, x, Linkage, Carriage)

n_beam = Linkage.ndof;
n_mass = 6;
n = n_beam + n_mass;

q  = x(1:n);
qd = x(n+1:2*n);

q_b    = q(1:n_beam);
q_mass = q(n_beam+1:end);

qd_b    = qd(1:n_beam);
qd_mass = qd(n_beam+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Solve geometric projection (ALGEBRAIC)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = ProjectS(Linkage, q_b, q_mass);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Kinematics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_s = FwdKinematicsAtS(Linkage, q_b, s);
T_m = variable_expmap_g(q_mass);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Beam + carriage dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ID, tau, ~, ~, ~, Mb] = DAEJacobians(Linkage, 0, q_b, qd_b, zeros(n_beam,1), [],[],[]);

Mm = body_inertia(Carriage.mass, ...
                   diag(Carriage.I), ...
                   Carriage.com_offset);

Cm = dinamico_coadj(qd_mass) * Mm;

gm = -Mm * dinamico_Adjoint(ginv(T_s)) * Linkage.G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Constraint residual + Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[err_mass, J_mass] = ErrorJAtS(Linkage, Carriage, s, q, q_mass);

% full constraint Jacobian
J = J_mass(:,1:n);

c = err_mass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Time derivative of constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Jdot_qdot = (J * qd - J * qd) / 1e-8;  % placeholder stable version
% (replace with your analytical Jdot if available)

gamma = Jdot_qdot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Build KKT system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = blkdiag(Mb, Mm);

h = [
    ID;
    Cm*qd_mass + gm
];

A = [
    M   J';
    J   zeros(size(J,1))
];

rhs = [
     - h;
    -gamma
];

sol = A \ rhs;

qdd = sol(1:n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Return ODE state derivative
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dx = [
    qd;
    qdd
];

end