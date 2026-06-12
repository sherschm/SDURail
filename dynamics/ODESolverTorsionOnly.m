
function [tvec_out, x_out, xdot_out] = ODESolverTorsionOnly(Linkage,Carriage,x0,tf,dt)
    mass = Carriage.mass;
    com_offset = Carriage.com_offset;
    I_carriage=Carriage.I;

    n_beam  = Linkage.ndof;
    n_mass  = Carriage.ndof;
    n       = n_beam + n_mass;    % number of position DOFs
    
    t_baum = 0.5;
    alpha =  2/t_baum; %50.0;%0.0;%% damping  0.0;%0.0;%
    beta =   1/(t_baum^2) ; %100.0;%0.0;%% stiffness 
    
    T_full_init = FwdKinematics(Linkage, x0(1:n_beam));
    T_L_fixed = T_full_init(end-3:end, :);
    T_0_fixed = T_full_init(1:4,:);

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
        s = ProjectS(Linkage, q_b, q_mass);
        
        % ----------------------------
        % Beam dynamics
        % ----------------------------
        [ID_no_Mqdd, tau, e, dID_dq, Cb, Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);
        
        %[ID, tau, e, dID_dq, Cb, Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,qdd_b,[],[],[]);

        % ----------------------------
        % kinematics at s
        % ----------------------------
        T_here = FwdKinematicsAtS(Linkage, q_b, s);

        [err0, J0, Jd0, errL, JL, JdL] = ErrorDynamicsAt0andL(Linkage,Carriage, q_b, qd_b, T_0_fixed, T_L_fixed);
        %[err_mass,J_mass_full,Jdqd_mass_full] = ErrorJAtS(Linkage, Carriage, s, q, qd);
        [err_mass,J_mass_full,Jdqd_mass_full] = ErrorJAtS(Linkage, Carriage, s, q, qd);


        % ----------------------------
        % Carriage Dynamics
        % ---------------------------
        Jm = eye(6);%SE3RightJacobianFromPose(q_mass);
        I_body = diag([I_carriage I_carriage I_carriage]);
        Mm =  Jm'* body_inertia(mass, I_body, com_offset) * Jm ;
        Cm = Jm'* dinamico_coadj(qd_mass) * Mm* Jm  ;
        gm =  Jm'* Mm * dinamico_Adjoint(ginv(T_here)) * Linkage.G;

        
        %gm = dinamico_Adjoint(ginv(T_here)) *Mm* Linkage.G;
        
        % =================================================
        % Jacobians
        % =================================================

        
        % ----------------------------
        % full constraint matrix
        % ----------------------------
        J_c = [J_mass_full; % free to move along x
                J0;
               JL];

        Jd_c = [Jdqd_mass_full;
                Jd0; %; % free to move along x
               JdL];

        
         err_pos = [err_mass;
                    err0;
                    errL];
        
        err_vel = J_c * qd;
        
        M_full = zeros(n,n);
        M_full(1:n_beam,1:n_beam) = Mb;
        M_full(n_beam+1:end,n_beam+1:end) = Mm;
        
        C_full = zeros(n,1);
        C_full(n_beam+1:end) = Cm*qd_mass;
        
        g_full = zeros(n,1);
        g_full(1:n_beam) = -ID_no_Mqdd + tau;
        g_full(n_beam+1:end) = gm;
        
        K_full = zeros(n,1);
        
        f = g_full - C_full - K_full;
        
        %gamma = -Jd_c*qd - alpha*err_vel - beta*err_pos;
        gamma =  - alpha*err_vel - beta*err_pos;

        nc = size(J_c,1);
        
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
        s = ProjectS(Linkage, q_b, q_mass);
        fprintf('s = %.3f\r', s);
    end
    
end