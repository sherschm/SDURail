
function [tvec_out, x_out, xdot_out] = BaumgarteBeamOnly(Linkage,x0,tf,dt)
    ndof = Linkage.ndof;

    n  = ndof;

    alpha = 0.0;% 0.0;%  % damping
    beta =  0.0;%0.0;% % stiffness
    
    T_full_init = FwdKinematics(Linkage, x0(1:ndof));
    T_L_init = T_full_init(end-3:end, :);
    T_0_init = T_full_init(1:4,:);

    %support = "pin-pin"; 
    support = "fixed-fixed"; 

    function err = fk_L_func(Linkage, Xpos, T_L_init, n)

        q_b = Xpos(1:n);
        
        T_L_full = FwdKinematics(Linkage, q_b);
        T_L = T_L_full(end-3:end,:);
        
        if support == "pin-pin"
         %T_L = T_L(end-3:end,:);
         %err = T_L_init(1:3,4) - T_L(1:3,4);
         err_full = piecewise_logmap(T_L_init) - piecewise_logmap(T_L);
         err = err_full(4:6);

        elseif support == "fixed-fixed"
         %dg = T_L_init \ T_L;
         %err = piecewise_logmap(dg);
         err = piecewise_logmap(T_L_init) - piecewise_logmap(T_L);
        end
    end
    
    function err = fk_0_func(Xpos, T_0_init)
    
        T_0 =variable_expmap_g(Xpos(1:6));
        
        if support == "pin-pin"
            %err = T_0_init(1:3,4) - T_0(1:3,4); 
            err_full = piecewise_logmap(T_0_init) - piecewise_logmap(T_0);
            err = err_full(4:6);

        elseif support == "fixed-fixed"

        %if fixed-fixed
         %dg = T_0_init \ T_0;
         %err = piecewise_logmap(dg);
         err = piecewise_logmap(T_0_init) - piecewise_logmap(T_0);
        end
    end

    function dX = dynamics(t, X)

        disp(['t = ', num2str(t)]);
        
        % ----------------------------
        % unpack state
        % ----------------------------
        xpos = X(1:n);
        xvel = X(n+1:2*n);
        
        q_b    = xpos(1:n);
        
        qd_b    = xvel(1:n);
        
        % ----------------------------
        % Beam dynamics
        % ----------------------------

        [ID,tau,e,dID_dq,Cb,Mb] = DAEJacobians(Linkage,0.0,q_b,qd_b,zeros(Linkage.ndof,1),[],[],[]);

        %Mb = GeneralizedMassMatrix(Linkage, q_b);
        %Cb = GeneralizedCoriolisMatrix(Linkage, q_b, qd_b);
        %gb = GeneralizedExternalForce(Linkage, q_b);
        
        %Kb = findK(Linkage);
        %Db = findD(Linkage);
        
        M_full = zeros(n, n);
        M_full(1:n, 1:n) = Mb;
        
        C_full = zeros(n,1);
        %C_full(1:n) = (Cb + Db) * qd_b;

        g_full = zeros(n,1);
        %g_full(1:n) = gb;
        g_full(1:n) = -ID + tau;    

        K_full = zeros(n,1);
        %K_full(1:n) = Kb * q_b;
        
        % =================================================
        % constraint wrappers
        % =================================================
        
        fk_L_wrapper = @(Xpos) fk_L_func(Linkage, Xpos, T_L_init, n);
        
        fk_0_wrapper = @(Xpos) fk_0_func(Xpos, T_0_init);
        
        % =================================================
        % Jacobians (requires AD or finite diff)
        % =================================================
        
        J_full = Jacobian(Linkage,q_b);
        if support == "fixed-fixed"
            J_L_full = J_full(end-5:end,:);
            J_0_full = J_full(1:6,:);
        elseif support == "pin-pin"
            J_L_full = J_full(end-2:end,:);
            J_0_full = J_full(4:6,:);
        end
        

        % ----------------------------
        % full constraint matrix
        % ----------------------------
        J_c = [ %J_mass_full ;
               %J_L_full];
               J_0_full];
        
        err_pos = [% fk_mass_err;
                   %fk_L_wrapper(xpos)];
                   fk_0_wrapper(xpos)];
        
        err_vel = J_c * xvel;
        
        % ----------------------------
        % constraint solve
        % ----------------------------
        Y = M_full \ J_c';
        S = J_c * Y;
        
        % regularisation
        if cond(S) > 1e12
            S = S + 1e-9 * eye(size(S,1));
        end
        
        rhs_lambda = J_c * (M_full \ (g_full - C_full - K_full)) ...
                     + alpha * err_vel + beta * err_pos;
        
        lambda = S \ rhs_lambda;
        
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

    % Condition function: returns 0 when trigger happens
    function [value, isterminal, direction] = term_condition(t, x)
        
        ndof = Linkage.ndof;
        s = x(ndof + 1);
        
        L = Linkage.VLinks(1).L;
        
        % event trigger condition
        value = min(s - 0.1, (L - 0.1) - s);
        
        % stop integration when event hits zero
        isterminal = 1;
        
        % detect both increasing/decreasing crossings
        direction = 0;
    
    end

    tspan = [0 tf];

    % ----------------------------
    % ODE options with event function
    % ----------------------------
    opts = odeset('Events', @term_condition);
    
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