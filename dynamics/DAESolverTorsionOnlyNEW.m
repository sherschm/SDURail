
function [tvec_out, x_out, xdot_out, lam_out] = DAESolverTorsionOnly(Linkage,Carriage,x0,tf,dt)
    mass = Carriage.mass;
    com_offset = Carriage.com_offset;
    I_mass=Carriage.I;

    n_beam  = Linkage.ndof;
    n_mass  = Carriage.ndof;
    n       = n_beam + n_mass;    % number of position DOFs
    
    %number of carriage constraints (3 for pose, 2 for lateral position)
    if Carriage.fixed == true
        nc_carriage = 6;
    else
        nc_carriage = 5;
    end
    %nc = nc_carriage + 12; % both rail ends fixed (needed when simulating lateral deformation
    %nc = nc_carriage + 11; % both rail ends fixed (needed when simulating lateral deformation
    nc = nc_carriage + 7; % 6 constraints for rail floating frame, plus one more for x-rotation of other end

    t_baum = 0.1;
    alpha =    1/t_baum; %50.0;%0.0;%% damping
    beta =   2/(t_baum^2) ; %100.0;%0.0;%% stiffness
    
    T_full_init = FwdKinematics(Linkage, x0(1:n_beam));
    T_L_fixed = T_full_init(end-3:end, :);
    T_0_fixed = T_full_init(1:4,:);

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
        [ID,tau,e,dID_dq,Cb,Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,qdd_b,[],[],[]);
        
        % ----------------------------
        % kinematics at s
        % ----------------------------
        T_here = FwdKinematicsAtS(Linkage, q_b, s);
        
        % ----------------------------
        % Carriage Dynamics
        % ---------------------------
        I_body = diag([I_mass I_mass I_mass]);
        Mm = body_inertia(mass, I_body, com_offset);
        Cm = dinamico_coadj(qd_mass) * Mm;
        gm = -Mm * dinamico_Adjoint(ginv(T_here)) * Linkage.G;
        
        % =================================================
        % Jacobians (requires AD or finite diff)
        % =================================================
        
        J_full = Jacobian(Linkage,q_b);
        Jd_full = Jacobiandot_corrected(Linkage,q_b,qd_b);
        
       % T_full = FwdKinematics(Linkage,q_b);
        %T_tip = T_full(end-3:end,:);
        
        %err =  FKErrorAtL(Linkage, q_b, T_L_fixed, n_beam);
        [J_L_full,Jd_L_full] = ErrorJAtL(Linkage, q_b, qd_b, T_L_fixed);
        %J_L_full = horzcat([J_L_full(end-5:end-3,:);J_L_full(end-1:end,:)],zeros(5,n_mass));
        %Jd_L_full = [Jd_L_full(end-5:end-3,:);Jd_L_full(end-1:end,:)];
        %J_L_full = horzcat(J_L_full(end-5:end,:),zeros(6,n_mass));
        %Jd_L_full = Jd_L_full(end-5:end,:);
        Jd_L_full = Jd_L_full(end-5,:);
        J_L_full = horzcat(J_L_full(end-5,:),zeros(1,n_mass));


        J_L_full = horzcat(J_full(end-5,:),zeros(1,n_mass));
        J_0_full = horzcat(J_full(1:6,:),zeros(6,n_mass));
        
        %Jd_L_full = horzcat(Jd_full(end-5:end,:),zeros(6,n_mass));
        Jd_L_full = horzcat(Jd_full(end-5,:),zeros(1,n_mass));
        Jd_0_full = horzcat(Jd_full(1:6,:),zeros(6,n_mass));
           
        % ----------------------------
        % mass constraint Jacobian
        % ----------------------------
        [J_mass_full, Jd_mass_full] = ErrorJAtS(Linkage, Carriage, s, q, qd);
        %J_mass_full = horzcat(-JacobianAtS(Linkage, q_b,s), JacobianAtPoint(Linkage, q_mass)); %fixes mass to point
        %Jd_mass_full = horzcat(-JacobianDotAtS(Linkage, q_b, qd_b,s), JacobiandotAtPoint(Linkage, q_mass,qd_mass)); %fixes mass to point
        
        if Carriage.fixed == true
            J_mass_corrected = J_mass_full(1:6,:);
            Jd_mass_corrected = Jd_mass_full(1:6,:);
        else
            %J_mass_corrected = [J_mass_full(1:3,:) ;J_mass_full(5:6,:)];
            J_mass_corrected = [J_mass_full(1:3,:) ;J_mass_full(5:6,:)]; %+ dphi_ds * ds_dq'; %correction
            Jd_mass_corrected = [Jd_mass_full(1:3,:) ;Jd_mass_full(5:6,:)]; % + dphi_ds * ds_dq'; %correction
        end

        % ----------------------------
        % full constraint matrix
        % ----------------------------
        FKErrorAtL_full = FKErrorAtL(Linkage, q_b, T_L_fixed, n_beam);
        %disp(FKErrorAtL_full(1))

        err_pos = [FKErrorAtS(Linkage, Carriage, s, q);
                    FKErrorAtL_full;
                   FKErrorAt0(q_b, T_0_fixed)];

        J_c = [J_mass_corrected; % free to move along x
               J_L_full;
               J_0_full];
        
        err_vel = J_c * q;

        Jd_cqd = [Jd_mass_corrected; % free to move along x
               Jd_L_full*qd;
               Jd_0_full*qd];

        
        F1 = q_dot - qd;
        
        %Rail Dynamics!
        F2 = ID - tau + J_c(:,1:n_beam)' * lambda;

        %carriage dynamics! (s
        F3 =  Mm*qdd_mass + Cm * qd_mass + gm + J_c(:,n_beam+1:n_beam+n_mass)' * lambda;% +[10000*t;0*t;0*t;0.0;0;0];
        
        %F4 = err_pos;
        %F4 = J_c * qd;
        %F4 = J_c * qddot + Jd_c*qd+ [Jd_qd; zeros(nc-5,1)];
        %F4 = J_c * qddot + [Jd_qd; Jd_L_full*qd ;Jd_0_full*qd ];
        %F4 = J_c * qddot + [Jd_qd; Jd_L_full*qd ;Jd_0_full*qd ];

        F4 = J_c * qddot + Jd_cqd; %+ alpha*err_vel + beta*err_pos;
        
        %disp(err_pos)

        %F4 = J_c * qddot + Jd_c*qd + alpha*err_vel + beta*err_pos;
       
        % if cond(J_c) > 1e6
        %     warning('J_c near singular at t=%g', t);
        % end
        % 
        % if cond(M_full) > 1e6
        %     warning('M_full ill-conditioned at t=%g', t);
        % end

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