
<<<<<<< HEAD
function [tvec_out, x_out, xdot_out, lam_out] = DAESolverTorsionOnly(Linkage,Carriage,x0,tf,dt)
    mass = Carriage.mass;
    com_offset = Carriage.com_offset;
    I_carriage=Carriage.I;
=======
function [tvec_out, x_out, xdot_out, lam_out] = DAESolverTorsionOnly(Linkage,y0,ydot0,nc,tf,dt,Carriage)
    % Carriage parameters
    mass = Carriage.mass;
    com_offset = Carriage.com_offset;
    I_carriage = Carriage.I;
>>>>>>> ae7fd69e1a59f6ad5da0ea14aa0656d761382ebe

    n_beam  = Linkage.ndof;
    n_mass  = Carriage.ndof;
    n       = n_beam + n_mass;    % number of position DOFs
    
<<<<<<< HEAD

    %number of carriage constraints (3 for pose, 2 for lateral position)
    if Carriage.fixed == true
        nc_carriage = 6;
    else
        nc_carriage = 5;
    end
    %nc_carriage = 6;

    %nc = nc_carriage + 12; % both rail ends fixed (needed when simulating lateral deformation
    %nc = nc_carriage + 11; % both rail ends fixed (needed when simulating lateral deformation
    nc = nc_carriage + 7; % 6 constraints for rail floating frame, plus one more for x-rotation of other end
    %nc = 12; % both rail ends fixed (needed when simulating lateral deformation

    t_baum = 0.1;
    alpha =    1/t_baum; %50.0;%% damping
    beta =   2/(t_baum^2) ; %100.0;%0.0;%% stiffness
    
    T_full_init = FwdKinematics(Linkage, x0(1:n_beam));
    T_L_fixed = T_full_init(end-3:end, :);
    T_0_fixed = T_full_init(1:4,:);

=======
    %Baumgarte stabilisation parameters 
    t_baum = 0.1; %Baumgarte time-constants
    alpha =    1/t_baum; %50.0;%0.0;%% damping
    beta =   2/(t_baum^2) ; %100.0;%0.0;%% stiffness
    
    T_full_init = FwdKinematics(Linkage, y0(1:n_beam));
    T_L_init = T_full_init(end-3:end, :);
    T_0_init = T_full_init(1:4,:);

    function err = carriage_constraint_err(Linkage, s, q)
        % Split generalized coordinates
        n_beam = Linkage.ndof;
        n_mass = 6;
        
        q_b = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        
        err_full = piecewise_logmap(FwdKinematicsAtS(Linkage,q_b,s))-q_mass;
        if Carriage.fixed == true
            err = err_full;
        else
            err = [err_full(1:3);err_full(5:6)];
        end
    end

   function err = fk_0_func(Xpos, T_0_init)
        %Function that describes error between floating frame coordinates
        %and fixed frame
        T_0 =variable_expmap_g(Xpos(1:6));
        err = piecewise_logmap(T_0_init) - piecewise_logmap(T_0);
    end

    function err = fk_L_func(Linkage, Xpos, T_L_init, n_beam)
        %Function that describes error between Rail end pose
        %and fixed pose
        q_b = Xpos(1:n_beam);
        
        T_L_full = FwdKinematics(Linkage, q_b);
        T_L = T_L_full(end-3:end,:);
        
         %dg = T_L_init \ T_L;
         %err = piecewise_logmap(dg);
         err_full = piecewise_logmap(T_L_init) - piecewise_logmap(T_L);
         err = err_full(1);

    end
   
>>>>>>> ae7fd69e1a59f6ad5da0ea14aa0656d761382ebe
    function F = dae_fun(t, y, ydot)
        q     = y(1:n);
        qd    = y(n+1:2*n);
        lambda = y(2*n+1:end); %Lagrange Multipliers
        
        qddot = ydot(n+1:2*n);
        q_dot  = ydot(1:n);

        % Parse State vector
        q_b    = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        
        qd_b    = qd(1:n_beam);
        qd_mass = qd(n_beam+1:n_beam+n_mass);
        
        qdd_b    = qddot(1:n_beam);
        qdd_mass = qddot(n_beam+1:n_beam+n_mass);

        % ----------------------------
        % solve internal coordinate s
        % ----------------------------
        s = ProjectS(Linkage, q_b, q_mass);
        
        % ----------------------------
        % Beam dynamics
        % ----------------------------
        %[ID_no_Mqdd, tau, e, dID_dq, Cb, Mb] = ...DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);

        [ID, tau, e, dID_dq, Cb, Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,qdd_b,[],[],[]);

        % ----------------------------
        % kinematics at s
        % ----------------------------
        T_here = FwdKinematicsAtS(Linkage, q_b, s);
        
        % ----------------------------
        % Carriage Dynamics
        % ---------------------------
        I_body = diag([I_carriage I_carriage I_carriage]);
        Mm = body_inertia(mass, I_body, com_offset);
        Cm = dinamico_coadj(qd_mass) * Mm;
        gm = -Mm * dinamico_Adjoint(ginv(T_here)) * Linkage.G;
        
        % =================================================
        % Jacobians
        % =================================================
        [err0, J0, Jd0, errL, JL, JdL] = ErrorDynamicsAt0andL(Linkage,Carriage, q_b, qd_b, T_0_fixed, T_L_fixed);

        [err_mass,J_mass_full,Jdqd_mass_full] = ErrorJAtS(Linkage, Carriage, s, q, qd);
        
        err_mass = [err_mass(1:3,:);err_mass(5:6,:)];
        J_mass_full = [J_mass_full(1:3,:);J_mass_full(5:6,:)];
        Jdqd_mass_full = [Jdqd_mass_full(1:3,:);Jdqd_mass_full(5:6,:)];
        
        
        % ----------------------------
        % full constraint matrix
        % ----------------------------
        J_c = [J_mass_full; % free to move along x
                J0;
               JL];

        Jd_c = [Jdqd_mass_full*qd;
                Jd0*qd; % free to move along x
               JdL*qd];
        
         err_pos = [err_mass;
                    err0;
                    errL
                   ];
        
        err_vel = J_c * qd;

        F1 = q_dot - qd;

        %F2 = M_full * qddot - (g_full - C_full - K_full - J_c' * lambda) +[zeros(Linkage.ndof,1);1000*t;0*t;0*t;0.0;0;0];
        
        %Rail Dynamics!
        %F2 = Mb*qdd_b + ID_no_Mqdd - tau + J_c(:,1:n_beam)' * lambda;
        F2 = ID - tau + J_c(:,1:n_beam)' * lambda;

        %carriage dynamics!
        F3 =  Mm*qdd_mass + Cm * qd_mass +gm + J_c(:,n_beam+1:n_beam+n_mass)' * lambda;
        
        %
        %F4 = err_pos;
        %F4 = err_vel;
        %F4 = J_c * qddot + Jd_c*qd+ [Jd_qd; zeros(nc-5,1)];
        %F4 = J_c * qddot + [Jd_qd; Jd_L_full*qd ;Jd_0_full*qd ];
        %F4 = J_c * qddot + [Jd_qd; Jd_L_full*qd ;Jd_0_full*qd ];

        %F4 = J_c * qddot + Jd_c*qd ;%+ alpha*err_vel + beta*err_pos;

        %F4 = J_c * qddot + Jd_c*qd + alpha*err_vel + beta*err_pos;
        F4 = J_c * qddot + Jd_c;%+ alpha*err_vel + beta*err_pos;


        F = [F1;
            F2;
            F3;
            F4];
    end

    lambda0_guess = zeros(nc,1);
    y0 = [x0;lambda0_guess];
    ydot0 = [x0(n+1:2*n);zeros(size(lambda0_guess,1)+n,1)];
      [y0_consistent, ydot0_consistent] = decic( ...
             @dae_fun, 0, ...
             y0, [true(n_beam,1); true(n_mass,1); true(n_beam,1); true(n_mass,1); false(nc,1)], ...
             ydot0, [false(n,1);false(n,1); true(nc,1)]);

    % 
    % ----------------------------
    % ODE options with event function
    % ----------------------------
    opts = odeset('OutputFcn',  @(t,y,flag) showTime(t,y,flag, Linkage));

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