function plotResults(tvec_out, x_sol, xd_sol, dt, Linkage, Carriage, T_0_fixed, T_L_fixed)

    % ----------------------------
    % Split state
    % ----------------------------
    q_beam_full = x_sol(:, 1:Linkage.ndof);
    q_mass_full = x_sol(:, Linkage.ndof+1:Linkage.ndof+6);

    N = length(tvec_out);

    % ----------------------------
    % Preallocate
    % ----------------------------
    s      = zeros(N,1);
    x_mass = zeros(N,1);
    y_mass = zeros(N,1);
    z_mass = zeros(N,1);

    % ----------------------------
    % Compute s and mass position
    % ----------------------------
    for i = 1:N
        qb = q_beam_full(i,:)';
        qm = q_mass_full(i,:)';

        % Project s
        s(i) = ProjectS(Linkage, qb, qm, Linkage.VLinks(1).L/2);  
        % (replace with your equivalent of project_s_newton if different)

        % Compute transform
        T = variable_expmap_g(qm);

        x_mass(i) = T(1,4);
        y_mass(i) = T(2,4);
        z_mass(i) = T(3,4);
    end

    % ----------------------------
    % Plot mass position
    % ----------------------------
    figure;
    plot(tvec_out, x_mass, 'DisplayName','x (m)', 'LineWidth', 2); hold on;
    plot(tvec_out, y_mass, 'DisplayName','y (m)', 'LineWidth', 2);
    plot(tvec_out, z_mass, 'DisplayName','z (m)', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Position (m)');
    legend;
    legend('Location','east');
    grid on;

    % Save
    if ~exist('plots','dir')
        mkdir('plots');
    end
    set(gcf, 'Units','pixels', 'Position',[100 100 500 350]);
    set(findall(gcf,'type','axes'),'FontSize',14);
    
    set(gcf, 'PaperPositionMode','auto');
    print('plots/response.png','-dpng','-r300');

    saveas(gcf, fullfile('plots','mass_response.png'));

    % ----------------------------
    % Plot s derivative
    % ----------------------------
    dsdt = diff(s) / dt;

    figure;
    plot(tvec_out(1:end-1), dsdt);
    xlabel('Time (s)');
    ylabel('ds/dt');
    grid on;

    % Save
    set(gcf, 'Units','pixels', 'Position',[100 100 500 350]);
    set(findall(gcf,'type','axes'),'FontSize',14);

    set(gcf, 'PaperPositionMode','auto');
    print('plots/sdot.png','-dpng','-r300');

    % ----------------------------
    % Combined plots
    % ----------------------------
    figure;

    % Left: s(t)
    subplot(1,2,1);
    plot(tvec_out, s, 'k');
    ylabel('s (m)');
    xlabel('Time (s)');
    grid on;

    % Right: beam coordinates
    subplot(1,2,2);
    plot(tvec_out, q_beam_full);
    ylabel('Beam modal strain coordinates');
    xlabel('Time (s)');
    grid on;

    % Save
    set(gcf, 'Units','pixels', 'Position',[100 100 1000 350]);
    set(findall(gcf,'type','axes'),'FontSize',14);
    
    set(gcf, 'PaperPositionMode','auto');
    print('plots/response.png','-dpng','-r300');

    % plot mass Constraint Violations
    ndof_tot = Linkage.ndof+Carriage.ndof;
    % Evaluate once to determine size of e_out
    %disp(x_sol(1, 1:ndof_tot))
    [err0,~,~,errL,~,~] = ErrorDynamicsAt0andL( ...
                 Linkage, Carriage, x_sol(1, 1:Linkage.ndof), zeros(Linkage.ndof,1), T_0_fixed, T_L_fixed);
    e_out = [err0;errL];
    n_e = length(e_out);
    
    E_out_all = zeros(n_e, N);
    
    % Loop over trajectory
    for k = 1:N
        [err0,~,~,errL,~,~] = ErrorDynamicsAt0andL( ...
                 Linkage, Carriage, x_sol(k, 1:Linkage.ndof), xd_sol(k, 1:Linkage.ndof), T_0_fixed, T_L_fixed);
        e_out = [err0;errL];
        E_out_all(:,k) = e_out(:);
    end
    
    % Time vector (assume unit timestep unless you have dt)
    t = 1:N;
    
    % Plot each constraint violation
    figure; hold on; grid on;
    colors = lines(n_e);
    
    for i = 1:n_e
        plot(t, E_out_all(i,:), 'LineWidth', 1.5, 'Color', colors(i,:));
    end
    
    xlabel('Time step');
    ylabel('Constraint violation (e\_out)');
    title('Constraint Violations Over Time');
    legend(arrayfun(@(i) sprintf('e_{out,%d}', i), 1:n_e, 'UniformOutput', false));
    print('plots/constraint_violation_ends.png','-dpng','-r300');
    

    % plot Constraint Violations
    ndof_tot = Linkage.ndof+Carriage.ndof;
    % Evaluate once to determine size of e_out

    [e_out, ~, ~] = ErrorJAtS(Linkage, Carriage, s(1), x_sol(1, 1:ndof_tot),xd_sol(1, 1:ndof_tot));
    n_e = length(e_out);
    
    E_out_all = zeros(n_e, N);
    
    % Loop over trajectory
    for k = 1:N
        [e_out, ~, ~] = ErrorJAtS(Linkage, Carriage, s(k), x_sol(k, 1:ndof_tot)',xd_sol(k, 1:ndof_tot)');
        E_out_all(:,k) = e_out(:);
    end
    
    % Time vector (assume unit timestep unless you have dt)
    t = 1:N;
    
    % Plot each constraint violation
    figure; hold on; grid on;
    colors = lines(n_e);
    
    for i = 1:n_e
        plot(t, E_out_all(i,:), 'LineWidth', 1.5, 'Color', colors(i,:));
    end
    
    xlabel('Time step');
    ylabel('Constraint violation (e\_out)');
    title('Constraint Violations Over Time');
    legend(arrayfun(@(i) sprintf('e_{out,%d}', i), 1:n_e, 'UniformOutput', false));
    print('plots/constraint_violation_mass.png','-dpng','-r300');
    
end