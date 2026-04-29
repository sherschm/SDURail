
function [tvec_out, x_out, xdot_out, lam_out] = DAESolverWithSVariable(Linkage,y0,ydot0,nc,tf,dt,support,fixed_carriage)
    ndof = Linkage.ndof;
    ndof_mass = 6;
    mass = 1.0;
    I_mass=0.0001;

    n_beam  = ndof;
    n_mass  = ndof_mass;
    n       = n_beam + n_mass + 1;    % number of position DOFs
    
    t_baum = 0.2;
    alpha =    1/t_baum; %50.0;%0.0;%% damping
    beta =   2/(t_baum^2) ; %100.0;%0.0;%% stiffness
    
    T_full_init = FwdKinematics(Linkage, y0(1:ndof));
    T_L_init = T_full_init(end-3:end, :);
    T_0_init = T_full_init(1:4,:);

    function err = carriage_constraint_err(Linkage, q)
        % Split generalized coordinates
        n_beam = Linkage.ndof;
        n_mass = 6;
        
        q_b = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        s = q(n);
        
        err_full = piecewise_logmap(FwdKinematicsAtS(Linkage,q_b,s))-q_mass;
        if fixed_carriage == true
            err = err_full;
        else
            err = [err_full(1:3);err_full(5:6)];
        end
    end

    function err = fk_L_func(Linkage, Xpos, T_L_init, n_beam)

        q_b = Xpos(1:n_beam);
        
        T_L_full = FwdKinematics(Linkage, q_b);
        T_L = T_L_full(end-3:end,:);
        
        if support == "pin-pin"
         T_L = T_L(end-3:end,:);
         err = T_L(1:3,4) - T_L_init(1:3,4);
         %err_full = piecewise_logmap(T_L_init) - piecewise_logmap(T_L);
         %err = err_full(4:6);

        elseif support == "fixed-fixed"
         %dg = T_L_init \ T_L;
         %err = piecewise_logmap(dg);
         err = piecewise_logmap(T_L_init) - piecewise_logmap(T_L);

        elseif support == "pin-fixed"
         %dg = T_L_init \ T_L;
         %err = piecewise_logmap(dg);
         err = piecewise_logmap(T_L_init) - piecewise_logmap(T_L);
        end
    end
    
    function err = fk_0_func(Xpos, T_0_init)
    
        T_0 =variable_expmap_g(Xpos(1:6));
        
        if support == "pin-pin"
            %err =  T_0(1:3,4)- T_0_init(1:3,4); 
            err_full = piecewise_logmap(T_0_init) - piecewise_logmap(T_0);
            err = err_full(4:6);

        elseif support == "fixed-fixed"

        %if fixed-fixed
         %dg = T_0_init \ T_0;
         %err = piecewise_logmap(dg);
         err = piecewise_logmap(T_0_init) - piecewise_logmap(T_0);

        elseif support == "pin-fixed"
            %err =  T_0(1:3,4) - T_0_init(1:3,4); 
            err_full = piecewise_logmap(T_0_init) - piecewise_logmap(T_0);
            err = err_full(4:6);
        end
    end
    

    function F = dae_fun(t, y, ydot)
        q     = y(1:n);
        qd    = y(n+1:2*n);
        lambda = y(2*n+1:2*n+nc);
        
        qddot = ydot(n+1:2*n);
        q_dot  = ydot(1:n);
        
        q_b    = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        s      = q(n);
        
        qd_b    = qd(1:n_beam);
        qd_mass = qd(n_beam+1:n_beam+n_mass);
        sd      = qd(n);
        
        % ----------------------------
        % Beam dynamics
        % ----------------------------

        [ID,tau,e,dID_dq,Cb,Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);
        
        % ----------------------------
        % kinematics at s
        % ----------------------------
        T_here = FwdKinematicsAtS(Linkage, q_b, s);
        
        % ----------------------------
        % rigid mass inertia
        % ---------------------------

        I_body = diag([I_mass I_mass I_mass]);
        Mm = body_inertia(mass, I_body, [0 0 0]);
        
        Cm = dinamico_coadj(qd_mass) * Mm;
        
        % ----------------------------
        % full system matrices
        % ----------------------------
        n_tot = n;
        
        M_full = zeros(n_tot, n_tot);
        M_full(1:n_beam, 1:n_beam) = Mb;
        M_full(n_beam+1:n_beam+n_mass, n_beam+1:n_beam+n_mass) = Mm;
        
        C_full = zeros(n_tot,1);
        %C_full(1:n_beam) = (Cb + Db) * qd_b;
        C_full(n_beam+1:n_beam+n_mass) = Cm * qd_mass;
        
        g_full = zeros(n_tot,1);
        %g_full(1:n_beam) = gb;
        g_full(1:n_beam) = -ID+tau;
        g_full(n_beam+1:n_beam+n_mass) = ...
            dinamico_Adjoint(ginv(T_here)) * mass * Linkage.G;
        
        K_full = zeros(n_tot,1);
        %K_full(1:n_beam) = Kb * q_b;
        
        % =================================================
        % constraint wrappers
        % =================================================
        
        fk_L_wrapper = @(Xpos) fk_L_func(Linkage, Xpos, T_L_init, n_beam);
        
        fk_0_wrapper = @(Xpos) fk_0_func(Xpos, T_0_init);
        
        % =================================================
        % Jacobians (requires AD or finite diff)
        % =================================================
        
        J_full = Jacobian(Linkage,q_b);
        Jd_full = Jacobiandot_corrected(Linkage,q_b,qd_b);
        
        T_full = FwdKinematics(Linkage,q_b);
        T_tip = T_full(end-3:end,:);

        J_end = dinamico_Adjoint(T_tip)*J_full(end-5:end,:);

        if support == "fixed-fixed"
            J_L_full = horzcat(J_full(end-5:end,:),zeros(6,n_mass+1));
            J_0_full = horzcat(J_full(1:6,:),zeros(6,n_mass+1));

            Jd_L_full = horzcat(Jd_full(end-5:end,:),zeros(6,n_mass+1));
            Jd_0_full = horzcat(Jd_full(1:6,:),zeros(6,n_mass+1));
        elseif support == "pin-pin"
            %J_L_full = horzcat(J_full(end-2:end,:),zeros(3,n_mass));
            J_L_full = horzcat(J_end(end-2:end,:),zeros(3,n_mass+1));
            J_0_full = horzcat(J_full(4:6,:),zeros(3,n_mass+1));

            Jd_L_full = horzcat(Jd_full(end-2:end,:),zeros(3,n_mass+1));
            Jd_0_full = horzcat(Jd_full(4:6,:),zeros(3,n_mass+1));
        elseif support == "pin-fixed"
            %J_L_full = horzcat(J_full(end-5:end,:),zeros(6,n_mass));
            J_L_full = horzcat(J_end(end-5:end,:),zeros(6,n_mass+1));
            J_0_full = horzcat(J_full(4:6,:),zeros(3,n_mass+1));
        end

        % ----------------------------
        % mass constraint Jacobian
        % ----------------------------

        J_mass_full = horzcat(JacobianAtS(Linkage, q_b,s),-eye(n_mass),zeros(n_mass,1)); %fixes mass to point
        Jd_mass_full = horzcat(JacobianDotAtS(Linkage, q_b, qd_b,s),zeros(n_mass,n_mass+1)); %fixes mass to point
        
        fk_mass_err = carriage_constraint_err(Linkage, q);
        %dfk_mass_ds = derivative_fd(@(ss) piecewise_logmap(FwdKinematicsAtS(Linkage,q_b,ss)), s);

        if fixed_carriage == true
            J_mass_corrected = J_mass_full(1:6,:);
            Jd_mass_corrected = Jd_mass_full(1:6,:);
        else
            %J_mass_corrected = [J_mass_full(1:3,:) ;J_mass_full(5:6,:)];
            J_mass_corrected = jacobian_fd(@(qq) carriage_constraint_err(Linkage, qq), q);
            Jd_mass_corrected = [Jd_mass_full(1:3,:) ;Jd_mass_full(5:6,:)];
        end

        % ----------------------------
        % full constraint matrix
        % ----------------------------
        J_c = [%J_mass_corrected; % free to move along x
               J_L_full;
               J_0_full];

        Jd_c = [Jd_mass_corrected; % free to move along x
               Jd_L_full;
               Jd_0_full];
        
        err_pos = [carriage_constraint_err(Linkage, q);% free to move along x
                   fk_L_wrapper(q);
                   fk_0_wrapper(q)];
        
        err_vel = J_c * qd;

        F1 = q_dot - qd;



        F2 = M_full * qddot - (g_full - C_full - K_full - J_c' * lambda);

        F3 = err_pos;
        %F3 = J_c * qd;
        %F3 = J_c * qddot + Jd_c*qd;
        %F3 = J_c * qddot + Jd_c*qd + alpha*err_vel + beta*err_pos;

        F = [F1;
             F2;
             F3];
        
    end

    tspan = [0 tf];
    

    [y0_consistent, ydot0_consistent] = decic( ...
    @dae_fun, 0, ...
    y0, [false(n,1); false(n,1); false(nc,1)], ...
    ydot0, [false(2*n,1); true(nc,1)]);
    
    % ----------------------------
    % ODE options with event function
    % ----------------------------
    %opts = odeset('RelTol',1e-4,'AbsTol',1e-6, 'OutputFcn',  @(t,y,flag) showTime(t,y,flag, Linkage));
    opts = odeset('OutputFcn',  @(t,y,flag) showTime(t,y,flag, Linkage));
    %opts = odeset('RelTol',1e-6,'AbsTol',1e-8, 'OutputFcn', @(t,y,flag) showTime(t,y,flag, Linkage));

    disp("Simulating...");

    [t,y] = ode15i(@dae_fun, [0 tf], y0_consistent, ydot0_consistent, opts);
   
    
    % ----------------------------
    % unpack solution
    % ----------------------------
    
    q_sol   = y(:,1:n);
    qd_sol  = y(:,n+1:2*n);
    lam_sol = y(:,2*n+1:end);
    
    tvec = t;
    
    % ----------------------------
    % interpolation to uniform grid
    % ----------------------------
    tvec_out = (0:dt:max(tvec))';
    
    x_out = interp1(tvec, q_sol, tvec_out, 'linear');
    xdot_out = interp1(tvec, qd_sol, tvec_out, 'linear');
    lam_out = interp1(tvec, lam_sol, tvec_out, 'linear');
end

function M = body_inertia(m, Inertia, com)

    % skew-symmetric matrix of COM
    px = [   0       -com(3)   com(2);
           com(3)      0      -com(1);
          -com(2)   com(1)      0    ];
    
    % composite inertia
    Ic = Inertia + m * (px * px');
    
    % spatial inertia matrix (6x6)
    M = [ Ic,      m * px;
          m * px', m * eye(3) ];
    
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