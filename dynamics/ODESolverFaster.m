
function [tvec_out, x_out, xdot_out] = ODEKKTSolver(Linkage,Carriage,x0,tf,dt,T_0_fixed,T_L_fixed, controller,carriage_constraints)

    n_beam  = Linkage.ndof;
    n_mass  = Carriage.ndof;
    n       = n_beam + n_mass;    % number of position DOFs
    
    %Baumgarte Stabilisation parameters

     t_baum = 0.1;
      alpha =  2/t_baum; %50.0;%0.0;%% damping 0.0;%0.0;%
      beta =   1/(t_baum) ; %100.0;%0.0;%% stiffness 
     % alpha = 0.0; %50.0;%0.0;%% damping 0.0;%0.0;%
     % beta =   0.0; %100.0;%0.0;%% stiffness 
     % W = blkdiag( ...
     %     10*eye(length(err_mass)), ...
     %     eye(length(err0)), ...
     %     eye(length(errL)));
    
     s_prev    = ProjectS(Linkage,Carriage, x0(1:n_beam), x0(n_beam+1:n_beam+n_mass),Linkage.VLinks(1).L/2);

    function xd = ode_fun(t, y)
        q     = y(1:n);
        qd    = y(n+1:2*n);
        
        % Parse State vector
        q_b    = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        
        qd_b    = qd(1:n_beam);
        qd_mass = qd(n_beam+1:n_beam+n_mass);

        % ----------------------------
        % solve internal coordinate s
        % ----------------------------
        tic
        [s,~] = ProjectS(...
            Linkage,...
            Carriage,...
            q_b,...
            q_mass,...
            s_prev);

        s_prev = s;

        t_proj = toc;
        %fprintf('t_proj = %e\n', t_proj);

        % ----------------------------
        % Beam dynamics
        % ----------------------------


        % kinematics
        %Compute constraint errors & constraint Jacobians
        tic
        [err0, J0, Jd0, errL, JL, JdL] = ErrorDynamicsAt0andL(Linkage,Carriage, q_b, qd_b, T_0_fixed, T_L_fixed);
        t_constraints= toc;
        %fprintf('t_end_constraints = %e\n', t_constraints);
        tic
        %
        if carriage_constraints == "both_ends"
            [err_mass, J_mass_full, Jd_mass_full] = ErrorJAtS_both_ends(Linkage, Carriage, s, q, qd);
        elseif carriage_constraints == "centre"
            [err_mass, J_mass_full, Jd_mass_full] = ErrorJAtS(Linkage, Carriage, s, q, qd);
        end
        % t_constraints= toc;
        %fprintf('t_s_constraints = %e\n', t_constraints);

        %Jd_mass_full = ErrorJdAtS_FD(Linkage, Carriage, s, q, qd) ; %Ideally replace with analytical Jdot.
        %Jd_mass_full = zeros(size(J_mass_full))

        % Carriage Dynamics
        %Mm =  GeneralizedMassMatrix(Carriage,q_mass);
        %Cm =  dinamico_coadj(qd_mass) * Mm  ;
        %gm =  Mm * dinamico_Adjoint(ginv(T_m)) * Linkage.G;
        % Build full constraint Jacobians
        J_c = [J_mass_full; % free to move along x
                J0;
               JL];

        Jd_c = [Jd_mass_full;
                Jd0; %; % free to move along x
               JdL];

        
        err_pos = [err_mass;
                    err0;
                    errL];
        
        err_vel = J_c * qd;
        
        %Build dynamics terms
        tic
        g_full = zeros(n,1);
        C_full = zeros(n,1);
        g_full = zeros(n,1);
        
        % [ID_no_Mqdd, tau, e, dID_dq, Cb, Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);
        % [ID_no_Mqddm, taum, em, dID_dqm, Cm, Mm] = DAEJacobians(Carriage,0.0,q_mass,qd_mass,zeros(Carriage.ndof,1),[],[],[]);
        % C_full(n_beam+1:end) = Cm*qd_mass;
        % M_full = zeros(n,n);
        % M_full(1:n_beam,1:n_beam) = Mb;
        % M_full(n_beam+1:end,n_beam+1:end) = Mm;
        % g_full(1:n_beam) = -ID_no_Mqdd + tau;
        % g_full(n_beam+1:end) = -ID_no_Mqddm + taum;


        M_full = zeros(n,n);
        M_full(1:n_beam,1:n_beam) = GeneralizedMassMatrix_corrected(Linkage,q_b);
        M_full(n_beam+1:end,n_beam+1:end) = GeneralizedMassMatrix_corrected(Carriage,q_mass);
        g_full(1:n_beam) = -dynamicsResidual( ...
            Linkage,...
            0.0,...
            [q_b;
            qd_b;
            zeros(n_beam,1)]);
        t_dyn= toc;
        %fprintf('t_dyn = %e\n', t_dyn);

        g_full(n_beam+1:n) = -dynamicsResidual( ...
                        Carriage,...
                        0.0,...
                        [q_mass;
                        qd_mass;
                        zeros(n_mass,1)]);
    

        % Control term
        % nu = 1;
        % Ix = [0;0;0;1;0;0];
        % Binb = [-JacobianAtS(Linkage,q_b,s)'*Ix ];
        % Binm = [Ix ];
        % B_full = [Binb;Binm];
        % PD controller!
        % example PD + geometry-based term
        %u = PD_carriage_control(Linkage, s, q, qd, controller);
        %u = 10000;
        %u=0.0;
    
        %Torsional force:

        % if t < 0.02
        %     torsion_tau = 1000000.0;
        % else
        %     torsion_tau = 0.0;
        % end
        % 
        % %Ix = [0;0;0;0;0;1];
        % Ix = [1;0;0;0;0;0];
        % Binb = zeros(n_beam,1);
        % Binm = Ix;
        % % Binb = JacobianAtS(Linkage,q_b,s)'*Ix;
        % % Binm = Ix;
        % 
        % B_full = [Binb;Binm];
        % 
        % f = g_full + B_full*torsion_tau;
        tic
        f = g_full;

        gamma = -Jd_c*qd - alpha*err_vel - beta*err_pos;
       
        % Construct the complete KKT matrix
        KKT = [M_full, J_c';
            J_c,    zeros(size(J_c,1))];

        rhs = [f; gamma];

        N = size(KKT,1);

        % Initialise diagonal scaling
        d = ones(N,1);

        % Ruiz equilibration
        maxIter = 10;

        for k = 1:maxIter

            D = spdiags(d,0,N,N);

            Ks = D * KKT * D;

            % Infinity norm of each row
            rowNorm = max(abs(Ks),[],2);

            % Avoid division by zero
            rowNorm = max(rowNorm,1e-12);

            % Update scaling factors
            alpha_ruiz = 1 ./ sqrt(rowNorm);

            % Limit individual scaling updates
            alpha_ruiz = max(min(alpha_ruiz,10),0.1);

            d = d .* alpha_ruiz;

        end

        % Final scaled system
        D = spdiags(d,0,N,N);

        KKTs = D * KKT * D;
        rhs_s = D * rhs;

        % fprintf('cond(KKT)  = %e\n', cond(KKT));
        % fprintf('cond(KKT2)  = %e\n', cond(KKTs));
        %fprintf('cond(Ms) = %e\n', cond(Ms));

        % Solve
        z = KKTs \ rhs_s;
        
        % Recover physical solution
        sol = D * z;

        qdd = sol(1:size(M_full,1));
        %lambda = sol(size(M_full,1)+1:end);

        xd = [qd;
              qdd];

        t_final= toc;
        %%fprintf('t_final = %e\n', t_final);
       
    end
    
    tspan = [0 tf];

    % ----------------------------
    % ODE options with event function
    % ----------------------------
    %opts = odeset('OutputFcn', @showTime,'Events', @term_condition);
    opts = odeset('OutputFcn',  @(t,y,flag) showTime(t,flag, s_prev));

    disp("Simulating...");
    
    [t, x] = ode15s(@ode_fun, tspan, x0, opts);
    %[t, x] = ode45(@ode_fun, tspan, x0, opts);

    
    % ----------------------------
    % unpack solution
    % ----------------------------
    
    x_sol    = x(:,1:n);
    xdot_sol = x(:,n+1:2*n);
    
    tvec = t;
    
    % ----------------------------
    % interpolation to uniform grid
    % ----------------------------
    tvec_out = (0:dt:max(tvec))';
    
    x_out = interp1(tvec, x_sol, tvec_out, 'linear');
    xdot_out = interp1(tvec, xdot_sol, tvec_out, 'linear');
end

function status = showTime(t,flag,s_prev)
    
    persistent next_print timer_start t_last
    
    status = 0;
    
    if strcmp(flag,'init')
    
        timer_start = tic;
        next_print = 1;
        t_last = [];
    
    elseif isempty(flag)
    
        elapsed = toc(timer_start);
    
        t_now = t(end);
    
        if isempty(t_last)
            h = NaN;
        else
            h = t_now - t_last;
        end
    
        while elapsed >= next_print
    
            fprintf(['Real time = %d s, ' ...
                'Sim time = %.6f s, ' ...
                'h = %.3e s, ' ...
                's = %.3f\n'], ...
                next_print, t_now, h, s_prev);
    
            next_print = next_print + 1;
        end
    
        t_last = t_now;
    
    elseif strcmp(flag,'done')
    
        clear t_last
    
    end

end
