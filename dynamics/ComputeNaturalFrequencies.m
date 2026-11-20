function f = ComputeNaturalFrequencies( ...
                    RailLinkage, RailCarriage, ...
                    s,q_guess, ...
                    T_0_fixed, T_L_fixed,carriage_constraints)

    n_beam = RailLinkage.ndof;
    n_mass = RailCarriage.ndof;
    n      = n_beam + n_mass;

    
    %compute equilibrium

    [q_equil, lambda0,~]  = FindEquilibriumFixed( ...
                    RailLinkage,RailCarriage, ...
                    s,q_guess, ...
                    T_0_fixed, T_L_fixed,carriage_constraints);

    disp(size(lambda0))

    q_b = q_equil(1:n_beam);
    q_mass = q_equil(n_beam+1:n);

    %Compute system dynamics
    g_full = zeros(n,1);
    M_full = zeros(n,n);
    [ID_no_Mqdd, tau, e, dID_dqb, Cb, Mb] = DAEJacobians(RailLinkage,0.0,q_b,zeros(n_beam,1),zeros(n_beam,1),[],[],[]);
    [ID_no_Mqddm, taum, em, dID_dqm, Cm, Mm] = DAEJacobians(RailCarriage,0.0,q_mass,zeros(n_mass,1),zeros(n_mass,1),[],[],[]);
    M_full = zeros(n,n);
    M_full(1:n_beam,1:n_beam) = Mb;
    M_full(n_beam+1:end,n_beam+1:end) = Mm;
    g_full(1:n_beam) = -ID_no_Mqdd + tau;
    g_full(n_beam+1:end) = -ID_no_Mqddm + taum;
    

    dID_dq = [dID_dqb zeros(n_beam, n_mass);
             zeros(n_mass,n_beam) dID_dqm];


    [err0, J0, Jd0, errL, JL, JdL] = ErrorDynamicsAt0andL(RailLinkage,RailCarriage, q_b, zeros(n_beam,1), T_0_fixed, T_L_fixed);

    if carriage_constraints == "both_ends"
        [err_mass, J_mass_full, Jd_mass_full] = ErrorJAtS_both_ends(RailLinkage, RailCarriage, s, q_equil, zeros(n,1));
    elseif carriage_constraints == "centre"
        [err_mass, J_mass_full, Jd_mass_full] = ErrorJAtS(RailLinkage, RailCarriage, s, q_equil, zeros(n,1));
    end
    
    J_c = [J_mass_full; J0; JL];
    % 
    % N = null(J_c,'r');
    % 
    % M_c = N'*M_full*N;
    % 
    % K_c = N'*dID_dq*N;   % check sign
    % 
    % [V,D] = eig(K_c,M_c);
    % 
    % lambda = real(diag(D));
    % 
    % tol = 1e-8;
    % lambda = lambda(lambda > tol);
    % 
    % f = sort(sqrt(lambda)/(2*pi));

    KKT = [M_full J_c';
       J_c    zeros(size(J_c,1))];

    rhs = [g_full;
           zeros(size(J_c,1),1)];
    
    sol = KKT \ rhs;
    
   % lambda0 = sol(n+1:end);

    h = 1e-7;

    K = zeros(n);
    
    for i = 1:n

        dq = zeros(n,1);
        dq(i) = h;
    
        qp = q_equil + dq;
        qm = q_equil - dq;
        
        [sp,~] = ProjectS(RailLinkage,RailCarriage,qp(1:n_beam),qp(n_beam+1:n),RailLinkage.VLinks(1).L/2);
        [sm,~] = ProjectS(RailLinkage,RailCarriage,qm(1:n_beam),qm(n_beam+1:n),RailLinkage.VLinks(1).L/2);

        Rp = Residual(sp, qp, lambda0,RailLinkage, RailCarriage, ...
                T_0_fixed, T_L_fixed, ...
                carriage_constraints);
        Rm = Residual(sm, qm, lambda0,RailLinkage, RailCarriage, ...
                T_0_fixed, T_L_fixed, ...
                carriage_constraints);
    
        K(:,i) = -(Rp - Rm)/(2*h);
    
    end
    disp(size(J_c))
    disp(norm(g_full + J_c'*lambda0))
    N  = null(J_c,'r');

    Mr = N'*M_full*N;
    Kr = N'*K*N;
    
    [V,D] = eig(Kr,Mr);
    
    omega = sqrt(real(diag(D)));
    f = sort(omega/(2*pi));
end

function R = Residual( ...
                s, q, lambda0, ...
                RailLinkage, RailCarriage, ...
                T_0_fixed, T_L_fixed, ...
                carriage_constraints)

    n_beam = RailLinkage.ndof;
    n_mass = RailCarriage.ndof;
    n      = n_beam + n_mass;

    q_b    = q(1:n_beam);
    q_mass = q(n_beam+1:end);

    qd_b    = zeros(n_beam,1);
    qd_mass = zeros(n_mass,1);

    %% Dynamics terms

    [ID_no_Mqdd, tau] = DAEJacobians( ...
        RailLinkage,0.0,q_b,qd_b,zeros(n_beam,1),[],[],[]);

    [ID_no_Mqddm, taum] = DAEJacobians( ...
        RailCarriage,0.0,q_mass,qd_mass,zeros(n_mass,1),[],[],[]);

    g_full = zeros(n,1);

    g_full(1:n_beam) = -ID_no_Mqdd  + tau;
    g_full(n_beam+1:end) = -ID_no_Mqddm + taum;

    %% Constraints

    [~, J0, ~, ~, JL, ~] = ...
        ErrorDynamicsAt0andL( ...
            RailLinkage,...
            RailCarriage,...
            q_b,...
            qd_b,...
            T_0_fixed,...
            T_L_fixed);

    if carriage_constraints == "both_ends"

        [~, J_mass_full, ~] = ...
            ErrorJAtS_both_ends( ...
                RailLinkage,...
                RailCarriage,...
                s,...
                q,...
                zeros(n,1));

    elseif carriage_constraints == "centre"

        [~, J_mass_full, ~] = ...
            ErrorJAtS( ...
                RailLinkage,...
                RailCarriage,...
                s,...
                q,...
                zeros(n,1));

    else

        error('Unknown carriage constraint type.');

    end

    J_c = [
        J_mass_full;
        J0;
        JL
    ];
    disp(size(J_c))
    %% Constrained residual

    R = g_full + J_c' * lambda0;

end