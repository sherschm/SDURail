
function [tvec_out, x_out, xdot_out] = Baumgarte(Linkage,x0,tf,dt,support)
    ndof = Linkage.ndof;
    ndof_mass = 6;
    mass = 1.0;
    I_mass=0.1;

    n_beam  = ndof;
    n_mass  = ndof_mass;
    n       = n_beam + n_mass;    % number of position DOFs
    
    t_baum = 0.1;
    alpha =    1/t_baum; %50.0;%0.0;%% damping
    beta =   2/(t_baum^2) ; %100.0;%0.0;%% stiffness
    
    T_full_init = FwdKinematics(Linkage, x0(1:ndof));
    T_L_init = T_full_init(end-3:end, :);
    T_0_init = T_full_init(1:4,:);

    function err = carriage_constraint_err(Linkage, s, q)
        % Split generalized coordinates
        n_beam = Linkage.ndof;
        n_mass = 6;
        
        q_b = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        
        err_full = piecewise_logmap(FwdKinematicsAtS(Linkage,q_b,s))-q_mass;
        err = [err_full(1:3);err_full(5:6)];
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

    function dX = dynamics(t, X)

        %disp(['t = ', num2str(t)]);
        
        % ----------------------------
        % unpack state
        % ----------------------------
        xpos = X(1:n);
        xvel = X(n+1:2*n);
        
        q_b    = xpos(1:n_beam);
        q_mass = xpos(n_beam+1:n_beam+n_mass);
        
        qd_b    = xvel(1:n_beam);
        qd_mass = xvel(n_beam+1:n_beam+n_mass);
        
        % ----------------------------
        % solve internal coordinate s
        % ----------------------------
        s = ProjectS(Linkage, q_b, q_mass);
        disp(['s = ', num2str(s)]);
        
        
        % ----------------------------
        % Beam dynamics
        % ----------------------------

        [ID,tau,e,dID_dq,Cb,Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);

        %Mb = GeneralizedMassMatrix(Linkage, q_b);
        %Cb = GeneralizedCoriolisMatrix(Linkage, q_b, qd_b);
        %gb = GeneralizedExternalForce(Linkage, q_b);
        
        %Kb = findK(Linkage);
        %Db = findD(Linkage);
        
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
        M_full(n_beam+1:end, n_beam+1:end) = Mm;
        
        C_full = zeros(n_tot,1);
        %C_full(1:n_beam) = (Cb + Db) * qd_b;
        C_full(n_beam+1:end) = Cm * qd_mass;
        
        g_full = zeros(n_tot,1);
        %g_full(1:n_beam) = gb;
        g_full(1:n_beam) = -ID+tau;
        g_full(n_beam+1:end) = ...
            dinamico_Adjoint(ginv(T_here)) * mass * Linkage.G;
        

        K_full = zeros(n_tot,1);
        %K_full(1:n_beam) = Kb * q_b;
        
        % =================================================
        % constraint wrappers
        % =================================================
        
        %fk_mass_wrapper = @(Xpos) carriage_constraint_err(Linkage, s, Xpos);
        
        fk_L_wrapper = @(Xpos) fk_L_func(Linkage, Xpos, T_L_init, n_beam);
        
        fk_0_wrapper = @(Xpos) fk_0_func(Xpos, T_0_init);
        
        % =================================================
        % Jacobians (requires AD or finite diff)
        % =================================================
        
        %J_L_full = jacobian_fd(fk_L_wrapper, xpos);
        %J_0_full = jacobian_fd(fk_0_wrapper, xpos);
        
        J_full = Jacobian(Linkage,q_b);
        Jd_full = Jacobiandot_corrected(Linkage,q_b,qd_b);
        
        T_full = FwdKinematics(Linkage,q_b);
        T_tip = T_full(end-3:end,:);

        J_end = dinamico_Adjoint(ginv(T_tip))*J_full(end-5:end,:);
        %Jd_end 
        %J_end = J_full(end-5:end,:);

        if support == "fixed-fixed"
            J_L_full = horzcat(J_full(end-5:end,:),zeros(6,n_mass));
            J_0_full = horzcat(J_full(1:6,:),zeros(6,n_mass));

            Jd_L_full = horzcat(Jd_full(end-5:end,:),zeros(6,n_mass));
            Jd_0_full = horzcat(Jd_full(1:6,:),zeros(6,n_mass));
        elseif support == "pin-pin"
            %J_L_full = horzcat(J_full(end-2:end,:),zeros(3,n_mass));
            J_L_full = horzcat(J_end(end-2:end,:),zeros(3,n_mass));
            J_0_full = horzcat(J_full(4:6,:),zeros(3,n_mass));

        elseif support == "pin-fixed"
            %J_L_full = horzcat(J_full(end-5:end,:),zeros(6,n_mass));
            J_L_full = horzcat(J_end(end-5:end,:),zeros(6,n_mass));
            J_0_full = horzcat(J_full(4:6,:),zeros(3,n_mass));
        end
        
        % ----------------------------
        % ds/dq
        % ----------------------------
        ds_dq = gradient_fd(@(q) ProjectS(Linkage, ...
                          q(1:n_beam), q(n_beam+1:end)), xpos);
        
        % ----------------------------
        % dφ/ds
        % ----------------------------
        dphi_ds = derivative_fd(@(ss) carriage_constraint_err(Linkage, ss, xpos), s);
        
        % ----------------------------
        % mass constraint Jacobian
        % ----------------------------
        %J_mass = jacobian_fd(fk_mass_wrapper, xpos);
        %J_mass_full = J_mass + dphi_ds * ds_dq';

        J_mass_full = horzcat(JacobianAtS(Linkage, q_b,s),-eye(n_mass)); %fixes mass to point
        
        Jd_mass_full = horzcat(JacobianDotAtS(Linkage, q_b, qd_b,s),zeros(n_mass,n_mass)); %fixes mass to point

        fk_mass_err = piecewise_logmap(FwdKinematicsAtS(Linkage,q_b,s))-q_mass;
        
        J_mass_corrected = [J_mass_full(1:3,:) ;J_mass_full(5:6,:)];% + dphi_ds * ds_dq'; %correction
        %J_mass_corrected = [J_mass_full(1:3,:) ;J_mass_full(5:6,:)] + dphi_ds * ds_dq'; %correction
        
        Jd_mass_corrected = [Jd_mass_full(1:3,:) ;Jd_mass_full(5:6,:)];% + dphi_ds * ds_dq'; %correction

        % ----------------------------
        % full constraint matrix
        % ----------------------------
        J_c = [J_mass_corrected % free to move along x
               J_L_full;
               J_0_full];
        
        Jd_c = [Jd_mass_corrected % free to move along x
               Jd_L_full;
               Jd_0_full];

        err_pos = [ fk_mass_err(1:3); fk_mass_err(5:6);% free to move along x
                   fk_L_wrapper(xpos);
                   fk_0_wrapper(xpos)];
        
        err_vel = J_c * xvel;
        
        % ----------------------------
        % constraint solve
        % ----------------------------
        %disp('M_full_cond = ');
        %display(cond(M_full));
        Y = M_full \ J_c';
        S = J_c * Y;
        
        % regularisation
        %if cond(S) > 1e12
        %    S = S + 1e-9 * eye(size(S,1));
        %end
        
        %rhs_lambda = J_c * (M_full \ (g_full - C_full - K_full)) ...
        %             + alpha * err_vel + beta * err_pos;
        rhs_lambda = J_c * (M_full \ (g_full - C_full - K_full)) ...
             + Jd_c * xvel ...
             + alpha * err_vel + beta * err_pos;
        
        %lambda = S \ rhs_lambda;
        lambda = pinv(S)*rhs_lambda;

        
        % ----------------------------
        % acceleration
        % ----------------------------
        qdd_full = M_full \ (g_full - C_full - K_full - J_c' * lambda);
        
        % ----------------------------
        % state derivative
        % ----------------------------
        dX = zeros(size(X));
        dX(1:n)    = xvel;
        dX(n+1:2*n) = qdd_full;
        
    end

    tspan = [0 tf];

    % ----------------------------
    % ODE options with event function
    % ----------------------------
    %opts = odeset('OutputFcn', @showTime,'Events', @term_condition);
    opts = odeset('OutputFcn', @showTime);

    disp("Simulating...");
    
    [t, x] = ode15s(@dynamics, tspan, x0, opts);
    
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

function status = showTime(t, y, flag)

    if isempty(flag)
        fprintf('t = %.3f\r', t(end));
    end

    status = 0;
end