
function [tvec_out, x_out, xdot_out] = ODEKKTSolverYOPMover(Linkage,Carriage,Mover,x0,tf,dt,T_0_fixed,T_L_fixed, controller)
    
    n_beam  = Linkage.ndof;
    n_mover  = 7;
    n       = n_beam + n_mover;    % number of position DOFs
    
    %Baumgarte Stabilisation parameters

    t_baum = 0.2;
    alpha =  2/t_baum; %50.0;%0.0;%% damping 0.0;%0.0;%
    beta =   1/(t_baum^2) ; %100.0;%0.0;%% stiffness 

    t_baum = 0.05;
    %alpha =  2/t_baum; %50.0;%0.0;%% damping 0.0;%0.0;%
    %beta =   1/(t_baum^2) ; %100.0;%0.0;%% stiffness 0.0;%
    
    alpha_terms =  [2*(2/t_baum)*ones(1,5) (2/t_baum)*ones(1,12)];
    beta_terms = [5*(1/(t_baum^2))*ones(1,5) (1/(t_baum^2))*ones(1,12)];

    alpha = diag(alpha_terms);
    beta = diag(beta_terms);

    function xd = ode_fun(t, y)
        q     = y(1:n);
        qd    = y(n+1:2*n);
        
        % Parse State vector
        q_b    = q(1:n_beam);
        q_mover = q(n_beam+1:n_beam+n_mover);
        
        qd_b    = qd(1:n_beam);
        qd_mover = qd(n_beam+1:n_beam+n_mover);

        % ----------------------------
        % solve internal coordinate s
        % ----------------------------
        [s,~] = ProjectS(Linkage, q_b, q_mover, Linkage.VLinks(1).L/2);
        
        % PD controller!
        % example PD + geometry-based term
        u = PD_carriage_control(Linkage, s, q, qd, controller);
        u_mover = [u;0;0;0;0;0];%
        % ----------------------------
        % Beam dynamics
        % ----------------------------
        [ID_no_Mqdd, tau, e, dID_dq, Cb, Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);

        %Mb = GeneralizedMassMatrix_corrected(Linkage,q_b);
        %Mb = GeneralizedMassMatrix(Linkage,q_b);

        % kinematics
        %T_rs = FwdKinematicsAtS(Linkage, q_b, s);
        T_m = variable_expmap_g(q_mover);

        %Compute constraint errors & constraint Jacobians
        [err0, J0, Jd0, errL, JL, JdL] = ErrorDynamicsAt0andL(Linkage, Carriage, q_b, qd_b, T_0_fixed, T_L_fixed);
        [err_mover, J_mover_full, Jd_mover_full] = ErrorJAtS_YOP(Linkage, Carriage, s, q, qd);

        %Jd_mover_full = ErrorJdAtS_FD(Linkage, Carriage, s, q, qd) ; %Ideally replace with analytical Jdot.
        %Jd_mover_full = zeros(size(J_mover_full))

        % Carriage Dynamics
        
        Mm_full =  assembleMassMatrix_mex(Mover.femodel.Mesh.Nodes, Mover.origins, Mover.Psi, Mover.DOF, Mover.m_i, Mover.weights, [q_mover;qd_mover]);
        Mm = Mm_full(8:end,8:end);

        fm_full = cb_dynamics_rhs(t, [q_mover;qd_mover], u_mover, Mover.rom.K, Mover.femodel.Mesh.Nodes, Mover.origins, Mover.Psi, Mover.DOF, Mover.m_i, Mover.weights);
        fm = fm_full(8:end);

        % Build full constraint Jacobians
        J_c = [J_mover_full; % free to move along x
                [J0 zeros(6,1)];
               [JL zeros(6,1)]];

        Jd_c = [Jd_mover_full;
                [Jd0 zeros(6,1)]; %; % free to move along x
               [JdL zeros(6,1)]];

        
        err_pos = [err_mover;
                    err0;
                    errL];
        
        err_vel = J_c * qd;
        
        %Build dynamics terms
       
        M_full = zeros(n,n);
        M_full(1:n_beam,1:n_beam) = Mb;
        M_full(n_beam+1:end,n_beam+1:end) = Mm;
        
        g_full = zeros(n,1);
        g_full(1:n_beam) = -ID_no_Mqdd + tau;
        g_full(n_beam+1:end) = fm;
        
        
        % Control term
        nu = 1;
        Ix = [0;0;0;1;0;0];
        Binb = [-JacobianAtS(Linkage,q_b,s)'*Ix ];
        Binm = [Ix ];
        B_full = [Binb;Binm];

        f = g_full + [Binb*u;u_mover;0];
        
        gamma = -Jd_c*qd - alpha*err_vel - beta*err_pos;

        nc = size(J_c,1);

        %Build KKT Matrix
        KKT = [M_full, J_c';
               J_c, zeros(nc)];
        %disp('KKT condition number:')
        %disp(cond(KKT))

        rhs = [f;
               gamma];
        
        sol = KKT \ rhs;
        
        qdd = sol(1:n);
        lambda = sol(n+1:end);
        
        xd = [qd;
              qdd];
       
    end
    
    tspan = [0 tf];

    % ----------------------------
    % ODE options with event function
    % ----------------------------
    %opts = odeset('OutputFcn', @showTime,'Events', @term_condition);
    opts = odeset('OutputFcn',  @(t,y,flag) showTime(t,y,flag, Linkage));

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

function status = showTime(t, y, flag, Linkage)

   
    %n_beam = 16; 
    n_beam = Linkage.ndof;

    n_mover = 6;
    %disp(['s = ', num2str(s)]);
    status = 0;
     if isempty(flag)
        fprintf('t = %.3f\r', t(end));
        q_b    = y(1:n_beam);
        q_mover = y(n_beam+1:n_beam+n_mover);

        % ----------------------------
        % solve internal coordinate s
        % ----------------------------
        [s,~] = ProjectS(Linkage, q_b, q_mover, Linkage.VLinks(1).L/2);
        fprintf('s = %.3f\r', s);
    end
    
end

% --- LOCAL FUNCTIONS ---
function rhs = cb_dynamics_rhs(t, Z, u, K, nodes, origins, Psi, DOF, m_i, weights)
    % IMPORTANT: This function now ONLY returns the right-hand side.
    
    F_in = [eye(6,6); zeros(1,6)];
    gbar = [zeros(7,6); F_in];
    
    % Just return f(Z) + g*u
    rhs = computeForce(K,nodes, origins, Psi, DOF, m_i, weights, Z) + gbar * u;
end

function Q_aug = computeForce(K,nodes, origin_ref, Psi, DOF, m_i, weights, Z)
    n = 7;
    Q_c= computeCoriolisForce_mex(nodes, origin_ref, Psi, DOF, m_i, weights, Z);
    Q_st = [zeros(n-1,1);-K(end,end) * 1e3 * Z(2*n)];
    %% M ddotx = Qc+Q_e+Qg+Qst
    Q = Q_c+Q_st;
    % 1. Extract angles and angular velocity
    euler_rates = Z(11:13);
    % 4. Assemble the corrected state derivative vector
    % [ translation_vel; euler_rates; modal_vel; Accelerations ]
    Q_aug = [Z(8:10); euler_rates; Z(14); Q];
end
