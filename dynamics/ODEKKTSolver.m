
function [tvec_out, x_out, xdot_out] = ODEKKTSolver(Linkage,Carriage,x0,tf,dt,T_0_fixed,T_L_fixed, controller)

    n_beam  = Linkage.ndof;
    n_mass  = Carriage.ndof;
    n       = n_beam + n_mass;    % number of position DOFs
    
    %Baumgarte Stabilisation parameters

     t_baum = 0.1;
     alpha =  2/t_baum; %50.0;%0.0;%% damping 0.0;%0.0;%
     beta =   1/(t_baum^2) ; %100.0;%0.0;%% stiffness 
     alpha = 0.0; %50.0;%0.0;%% damping 0.0;%0.0;%
     beta =   0.0; %100.0;%0.0;%% stiffness 
     % W = blkdiag( ...
     %     10*eye(length(err_mass)), ...
     %     eye(length(err0)), ...
     %     eye(length(errL)));


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
        [s,~] = ProjectS(Linkage,Carriage, q_b, q_mass, Linkage.VLinks(1).L/2);

        % ----------------------------
        % Beam dynamics
        % ----------------------------


        % kinematics
        %T_rs = FwdKinematicsAtS(Linkage, q_b, s);
        T_m = variable_expmap_g(q_mass);
        
        %Compute constraint errors & constraint Jacobians
        [err0, J0, Jd0, errL, JL, JdL] = ErrorDynamicsAt0andL(Linkage,Carriage, q_b, qd_b, T_0_fixed, T_L_fixed);
        % [err_mass, J_mass_full, Jd_mass_full] = ErrorJAtS(Linkage, Carriage, s, q, qd);
        [err_mass, J_mass_full, Jd_mass_full] = ErrorJAtS_both_ends(Linkage, Carriage, s, q, qd);

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

        g_full(n_beam+1:n) = -dynamicsResidual( ...
                        Carriage,...
                        0.0,...
                        [q_mass;
                        qd_mass;
                        zeros(n_mass,1)]);
    

        % Control term
        nu = 1;
        Ix = [0;0;0;1;0;0];
        Binb = [-JacobianAtS(Linkage,q_b,s)'*Ix ];
        Binm = [Ix ];
        B_full = [Binb;Binm];

        % PD controller!
        % example PD + geometry-based term
        %u = PD_carriage_control(Linkage, s, q, qd, controller);
        %u = 10000;
        u=0.0;

        f = g_full + B_full*u;
        
        gamma = -Jd_c*qd - alpha*err_vel - beta*err_pos;
       % gamma =  - alpha*err_vel - beta*err_pos;

        nc = size(J_c,1);
        
        % Scale coordinates so that diagonal(M) ~= 1
        d = diag(M_full);

        d = max(d, 1e-6*max(d));

        S = diag(1./sqrt(d));
        Ms = S' * M_full * S;

        %fprintf('cond(M)  = %e\n', cond(M_full));
        %fprintf('cond(Ms) = %e\n', cond(Ms));

        Js = J_c * S;

        KKTs = [Ms, Js';
            Js, zeros(size(Js,1))];
        %fprintf('cond(KKTs) = %e\n', cond(KKTs));
        rhs_s = [S' * f;
            gamma];

        sol_s = KKTs \ rhs_s;

        zdd    = sol_s(1:n);
        lambda = sol_s(n+1:end);

        qdd = S * zdd;

        xd = [qd;
              qdd];
       
    end
    
    tspan = [0 tf];

    % ----------------------------
    % ODE options with event function
    % ----------------------------
    %opts = odeset('OutputFcn', @showTime,'Events', @term_condition);
    opts = odeset('OutputFcn',  @(t,y,flag) showTime(t,y,flag, Linkage,Carriage));

    disp("Simulating...");
    
    [t, x] = ode15s(@ode_fun, tspan, x0, opts);
    
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

function status = showTime(t, y, flag, Linkage,Carriage)

   
    %n_beam = 16; 
    n_beam = Linkage.ndof;

    n_mass = 6;
    %disp(['s = ', num2str(s)]);
    status = 0;
     if isempty(flag)
        fprintf('t = %.3f\r', t(end));
        q_b    = y(1:n_beam);
        q_mass = y(n_beam+1:n_beam+n_mass);

        % ----------------------------
        % solve internal coordinate s
        % ----------------------------
        [s,~] = ProjectS(Linkage,Carriage, q_b, q_mass, Linkage.VLinks(1).L/2);
        fprintf('s = %.3f\r', s);
    end
    
end