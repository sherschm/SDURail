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
        
         err_pos = [err_mass;
                    err0;
                    errL
                   ];


        F1 = ID - tau + J_c(:,1:n_beam)' * lambda;

        %carriage dynamics!
        F2 =  Mm*qdd_mass + Cm * qd_mass +gm + J_c(:,n_beam+1:n_beam+n_mass)' * lambda;
        F3 = err_pos;
        
        %Solve with Ipopt
    end
