function [tvec_out,x_out,xdot_out] = ...
    NewmarkKKTSolver( ...
    Linkage,...
    Carriage,...
    x0,...
    tf,...
    dt,...
    T_0_fixed,...
    T_L_fixed,...
    controller,...
    carriage_constraints)

n_beam = Linkage.ndof;
n_mass = Carriage.ndof;
n      = n_beam + n_mass;

% Newmark parameters
gamma_nm = 0.5;
beta_nm  = 0.25;

Nt = floor(tf/dt) + 1;

tvec_out = (0:dt:tf)';

x_out    = zeros(Nt,n);
xdot_out = zeros(Nt,n);

% Initial conditions
q  = x0(1:n);
qd = x0(n+1:2*n);

s_prev = ProjectS( ...
    Linkage,...
    Carriage,...
    q(1:n_beam),...
    q(n_beam+1:end),...
    Linkage.VLinks(1).L/2);

% Initial acceleration
[qdd,s_prev] = ComputeAcceleration( ...
    q,...
    qd,...
    s_prev);

x_out(1,:)    = q.';
xdot_out(1,:) = qd.';

fprintf('Simulating...\n');

for k = 1:Nt-1

    t = tvec_out(k);

    %
    % Predictor
    %

    q_pred = q ...
        + dt*qd ...
        + dt^2*(0.5-beta_nm)*qdd;

    qd_pred = qd ...
        + dt*(1-gamma_nm)*qdd;

    %
    % Solve constrained acceleration
    %

    [qdd_new,s_prev] = ComputeAcceleration( ...
        q_pred,...
        qd_pred,...
        s_prev);

    %
    % Corrector
    %

    q_new = q_pred ...
        + beta_nm*dt^2*qdd_new;

    qd_new = qd_pred ...
        + gamma_nm*dt*qdd_new;

    %
    % Optional constraint projection
    %

    % [q_new,qd_new,s_prev] =
    %     ProjectConstraints(q_new,qd_new,s_prev);

    %
    % Store
    %

    q   = q_new;
    qd  = qd_new;
    qdd = qdd_new;

    x_out(k+1,:)    = q.';
    xdot_out(k+1,:) = qd.';

    if mod(k,100)==0
        fprintf('t = %.3f / %.3f\r',t,tf);
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qdd,s] = ComputeAcceleration(q,qd,s_guess)

    %
    % Parse states
    %

    q_b    = q(1:n_beam);
    q_mass = q(n_beam+1:end);

    qd_b    = qd(1:n_beam);
    qd_mass = qd(n_beam+1:end);

    %
    % Find carriage coordinate
    %

    [s,~] = ProjectS( ...
        Linkage,...
        Carriage,...
        q_b,...
        q_mass,...
        s_guess);

    %
    % Constraints
    %

    [err0,J0,Jd0,errL,JL,JdL] = ...
        ErrorDynamicsAt0andL( ...
        Linkage,...
        Carriage,...
        q_b,...
        qd_b,...
        T_0_fixed,...
        T_L_fixed);

    if carriage_constraints == "both_ends"

        [err_mass,J_mass_full,Jd_mass_full] = ...
            ErrorJAtS_both_ends( ...
            Linkage,...
            Carriage,...
            s,...
            q,...
            qd);

    else

        [err_mass,J_mass_full,Jd_mass_full] = ...
            ErrorJAtS( ...
            Linkage,...
            Carriage,...
            s,...
            q,...
            qd);

    end

    J_c = [
        J_mass_full
        J0
        JL
        ];

    Jd_c = [
        Jd_mass_full
        Jd0
        JdL
        ];

    %
    % Mass matrix
    %

    M_full = zeros(n,n);

    M_full(1:n_beam,1:n_beam) = ...
        GeneralizedMassMatrix_corrected( ...
        Linkage,q_b);

    M_full(n_beam+1:end,n_beam+1:end) = ...
        GeneralizedMassMatrix_corrected( ...
        Carriage,q_mass);

    %
    % Dynamics
    %

    g_full = zeros(n,1);

    g_full(1:n_beam) = ...
        -dynamicsResidual( ...
        Linkage,...
        0.0,...
        [q_b;
         qd_b;
         zeros(n_beam,1)]);

    g_full(n_beam+1:end) = ...
        -dynamicsResidual( ...
        Carriage,...
        0.0,...
        [q_mass;
         qd_mass;
         zeros(n_mass,1)]);

    f = g_full;

    %
    % Constraint acceleration equation
    %

    gamma = -Jd_c*qd;

    %
    % Coordinate scaling
    %

    d = diag(M_full);

    d = max(d,1e-6*max(d));

    S  = diag(1./sqrt(d));

    Ms = S'*M_full*S;

    Js = J_c*S;

    Minvf = Ms\(S'*f);

    A = Js*(Ms\Js');

    b = Js*Minvf - gamma;

    lambda = A\b;

    zdd = Minvf - Ms\(Js'*lambda);

    qdd = S*zdd;

end

end