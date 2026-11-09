function [qeq,qdeq,qddeq,lambdaeq,exitflag] = FindEquilibrium( ...
                    Linkage,Carriage, ...
                    s_guess,q_guess, ...
                    T_0_fixed,T_L_fixed)

    n_beam = Linkage.ndof;
    n_mass = Carriage.ndof;
    n      = n_beam + n_mass;

    %-------------------------------------------------------------
    % Determine number of constraints
    %-------------------------------------------------------------
    qd0 = zeros(n,1);

    [err0,~,~,errL,~,~] = ErrorDynamicsAt0andL( ...
                     Linkage,Carriage,...
                     q_guess(1:n_beam),...
                     zeros(n_beam,1),...
                     T_0_fixed,T_L_fixed);

    % [err_mass,~] = ErrorJAtS_fixed( ...
    %                  Linkage,Carriage,...
    %                  s,q_guess,qd0);
    [err_mass,~] = ErrorJAtS( ...
        Linkage,Carriage,...
        s_guess,q_guess,qd0);

    nc = length(err_mass) + length(err0) + length(errL);

    %-------------------------------------------------------------
    % Decision variables:
    %
    % x = [ q
    %       qd
    %       qdd
    %       lambda ]
    %-------------------------------------------------------------
    x0 = [q_guess;
          zeros(n,1);
          zeros(n,1);
          zeros(nc,1)];

    opts = optimoptions('fsolve',...
        'Algorithm','trust-region-dogleg',...
        'Display','iter',...
        'FunctionTolerance',1e-12,...
        'StepTolerance',1e-12,...
        'OptimalityTolerance',1e-12);

    [xsol,~,exitflag] = fsolve(@residual,x0,opts);

    qeq      = xsol(1:n);
    qdeq     = xsol(n+1:2*n);
    qddeq    = xsol(2*n+1:3*n);
    lambdaeq = xsol(3*n+1:end);

    %=============================================================
    % Residual
    %=============================================================
    function r = residual(x)
        
        q      = x(1:n);
        qd     = x(n+1:2*n);
        qdd    = x(2*n+1:3*n);
        lambda = x(3*n+1:end);

        q_b    = q(1:n_beam);
        qd_b   = qd(1:n_beam);
        qdd_b  = qdd(1:n_beam);

        q_mass   = q(n_beam+1:end);
        qd_mass  = qd(n_beam+1:end);
        qdd_mass = qdd(n_beam+1:end);

        [s,~] = ProjectS(Linkage,q_b,q_mass,s_guess);

        %---------------------------------------------------------
        % Constraint equations
        %---------------------------------------------------------
        [err0,J0,~, ...
         errL,JL,~] = ErrorDynamicsAt0andL( ...
                            Linkage,...
                            Carriage,...
                            q_b,...
                            qd_b,...
                            T_0_fixed,...
                            T_L_fixed);

        % [err_mass,...
        %  J_mass_full,...
        %  ~] = ErrorJAtS_fixed( ...
        %                     Linkage,...
        %                     Carriage,...
        %                     s,...
        %                     q,...
        %                     qd);
        [err_mass,...
            J_mass_full,...
            ~] = ErrorJAtS( ...
            Linkage,...
            Carriage,...
            s,...
            q,...
            qd);

        Jc = [J_mass_full;
              J0;
              JL];

        phi = [err_mass;
               err0;
               errL];

        %---------------------------------------------------------
        % Dynamics residuals
        %
        % g_full = M*qdd + n(q,qd) - Bu
        %---------------------------------------------------------
        g_full = zeros(n,1);

        g_full(1:n_beam) = dynamicsResidual( ...
                                Linkage,...
                                0.0,...
                                [q_b;
                                 qd_b;
                                 qdd_b]);

        g_full(n_beam+1:n) = dynamicsResidual( ...
                                Carriage,...
                                0.0,...
                                [q_mass;
                                 qd_mass;
                                 qdd_mass]);

        %---------------------------------------------------------
        % Force balance
        %
        % M*qdd + n - Bu + Jc'*lambda = 0
        %---------------------------------------------------------
        force_balance = g_full + Jc'*lambda;

        r = [force_balance;
            phi;
            qd;
            qdd];

    end

end
