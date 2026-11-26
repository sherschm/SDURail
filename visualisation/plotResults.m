function plotResults(tvec_out, x_out, xd_sol, dt, RailLinkage, RailCarriage, T_0_fixed, T_L_fixed,carriage_constraints)
    close all
    % ----------------------------
    % Split state
    % ----------------------------
    q_beam_full = x_out(:, 1:RailLinkage.ndof);
    q_mass_full = x_out(:, RailLinkage.ndof+1:RailLinkage.ndof+6);

    N = length(tvec_out);
    styles = {'-','--',':','-.'};   % add more if needed
    c_green = [0.0000 0.6000 0.3000];

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
        s(i) = ProjectS(RailLinkage,RailCarriage, qb, qm, RailLinkage.VLinks(1).L/2);  
        % (replace with your equivalent of project_s_newton if different)
        T = FwdKinematicsAtS(RailLinkage,qb, s(i));
        % Compute transform
        %T = variable_expmap_g(qm);

        x_mass(i) = T(1,4);
        y_mass(i) = T(2,4);
        z_mass(i) = T(3,4);
    end

    % ----------------------------
    % Plot mass position
    % ----------------------------
    figure('Color','w');

    t = tiledlayout(2,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');
    nexttile;

    plot(tvec_out, x_mass-x_mass(1),styles{2}, 'DisplayName','x', 'LineWidth', 2); 
    hold on;
    plot(tvec_out, y_mass-y_mass(1),styles{3}, 'DisplayName','y', 'LineWidth', 2);
    plot(tvec_out, z_mass-z_mass(1),styles{1}, 'DisplayName','z', 'LineWidth', 2,'Color', c_green);
    xlabel('Time (s)');
    ylabel('Carriage position (m)');
    legend;
    legend('Location','best', ...
        'Interpreter','latex', TextColor ='black');
    legend('Location','eastoutside');
    set(legend,'Color','w')
    set(legend,'TextColor','k')
    set(legend,'EdgeColor','k')
    grid on;
    title('(a)','Interpreter','latex','Interpreter','latex', ...
        'Color','k');
xlim([0,0.2])
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.FontSize = 11;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    grid on;
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.08;

    nexttile

    plot(tvec_out, q_mass_full(:,end-5),styles{2}, 'DisplayName','Roll', 'LineWidth', 2); 
    hold on;
    plot(tvec_out, q_mass_full(:,end-4),styles{3}, 'DisplayName','Pitch', 'LineWidth', 2);
    plot(tvec_out, q_mass_full(:,end-3),styles{1}, 'DisplayName','Yaw', 'LineWidth', 2,'Color', c_green);
    xlabel('Time (s)');
    ylabel('Carriage orientation (rad)');
    legend;
    legend('Location','best', ...
        'Interpreter','latex', TextColor ='black');
    legend('Location','eastoutside');
    set(legend,'Color','w')
    set(legend,'TextColor','k')
    set(legend,'EdgeColor','k')
    grid on;
    title('(b)','Interpreter','latex','Interpreter','latex', ...
        'Color','k');

    xlim([0,0.2])
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.FontSize = 11;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    grid on;
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.08;

    % Save
    if ~exist('plots','dir')
        mkdir('plots');
    end
    set(gcf,'Units','inches','Position',[1 1 10 5.5]);
    set(findall(gcf,'type','axes'),'FontSize',14);
    
    set(gcf, 'PaperPositionMode','auto');
    print('plots/response.png','-dpng','-r300');

    exportgraphics(gcf, ...
        'plots/carriage_response.pdf', ...
        'ContentType','vector');


    %% Forcing:

    % Time vector
    t = linspace(0,0.2,1000);

    % Torque pulse
    torsion_tau = zeros(size(t));
    torsion_tau(t < 0.02) = 10.0;

    % Plot
    figure('Color','w')
    plot(t,torsion_tau,'LineWidth',2,Color ='r')

    xlabel('Time (s)')
    ylabel('Torsional torque, $\tau$ (Nm)', ...
        'Interpreter','latex')
    title('Applied Torsional Pulse',Color = 'k')
    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.FontSize = 14;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.Box = 'off';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    grid on;
    ax.GridAlpha = 0.15;
    ax.MinorGridAlpha = 0.08;
    grid on
    box off
    set(gcf,'Units','inches','Position',[1 1 5 4]);

    xlim([0 0.2])
    ylim([0 11])

    exportgraphics(gcf, ...
        'plots/torsion_pulse.pdf', ...
        'ContentType','vector');
% %%
%     figure('Color','w');
% 
%     plot(tvec_out, z_mass, 'DisplayName','z (m)', 'LineWidth', 2);
%     xlabel('Time (s)');
%     ylabel('Position (m)');
%     saveas(gcf, fullfile('plots','mass_response_z_only.png'));
% 
%     % ----------------------------
%     % Plot s derivative
%     % ----------------------------
%     dsdt = diff(s) / dt;
% 
%     figure('Color','w');
%     plot(tvec_out(1:end-1), dsdt);
%     xlabel('Time (s)');
%     ylabel('ds/dt');
%     grid on;
% 
%     % Save
%     set(gcf, 'Units','pixels', 'Position',[100 100 500 350]);
%     set(findall(gcf,'type','axes'),'FontSize',14);
% 
%     set(gcf, 'PaperPositionMode','auto');
%     print('plots/sdot.png','-dpng','-r300');
% 
%     % ----------------------------
%     % Combined plots
%     % ----------------------------
%     figure('Color','w');
% 
%     % Left: s(t)
%     subplot(1,2,1);
%     plot(tvec_out, s, 'k');
%     ylabel('s (m)');
%     xlabel('Time (s)');
%     grid on;
% 
%     % Right: beam coordinates
%     subplot(1,2,2);
%     plot(tvec_out, q_beam_full);
%     ylabel('Beam modal strain coordinates');
%     xlabel('Time (s)');
%     grid on;
% 
%     % Save
%     set(gcf, 'Units','pixels', 'Position',[100 100 1000 350]);
%     set(findall(gcf,'type','axes'),'FontSize',14);
% 
%     set(gcf, 'PaperPositionMode','auto');
%     print('plots/response.png','-dpng','-r300');
% 
%     % plot mass Constraint Violations
%     ndof_tot = RailLinkage.ndof+RailCarriage.ndof;
%     % Evaluate once to determine size of e_out
%     %disp(x_out(1, 1:ndof_tot))
%     [err0,~,~,errL,~,~] = ErrorDynamicsAt0andL( ...
%                  RailLinkage, RailCarriage, x_out(1, 1:RailLinkage.ndof), zeros(RailLinkage.ndof,1), T_0_fixed, T_L_fixed);
%     e_out = [err0;errL];
%     n_e = length(e_out);
% 
%     E_out_all = zeros(n_e, N);
% 
%     % Loop over trajectory
%     for k = 1:N
%         [err0,~,~,errL,~,~] = ErrorDynamicsAt0andL( ...
%                  RailLinkage, RailCarriage, x_out(k, 1:RailLinkage.ndof), xd_sol(k, 1:RailLinkage.ndof), T_0_fixed, T_L_fixed);
%         e_out = [err0;errL];
%         E_out_all(:,k) = e_out(:);
%     end
% 
%     % Time vector (assume unit timestep unless you have dt)
%     t = 1:N;
% 
%     % Plot each constraint violation
%     figure('Color','w'); hold on; grid on;
%     colors = lines(n_e);
% 
%     for i = 1:n_e
%         plot(t, E_out_all(i,:), 'LineWidth', 1.5, 'Color', colors(i,:));
%     end
% 
%     xlabel('Time step');
%     ylabel('Constraint violation (e\_out)');
%     title('Constraint Violations Over Time');
%     legend(arrayfun(@(i) sprintf('e_{out,%d}', i), 1:n_e, 'UniformOutput', false));
%     print('plots/constraint_violation_ends.png','-dpng','-r300');
% 
% 
%     % plot Constraint Violations
%     ndof_tot = RailLinkage.ndof+RailCarriage.ndof;
%     % Evaluate once to determine size of e_out
% 
%     if carriage_constraints == "centre"
%         [e_out, ~, ~] = ErrorJAtS(RailLinkage, RailCarriage, s(1), x_out(1, 1:ndof_tot),xd_sol(1, 1:ndof_tot));
%     elseif carriage_constraints == "both_ends"
%         [e_out, ~, ~] = ErrorJAtS_both_ends(RailLinkage, RailCarriage, s(1), x_out(1, 1:ndof_tot)',xd_sol(1, 1:ndof_tot)');
%     end
%     n_e = length(e_out);
% 
%     E_out_all = zeros(n_e, N);
% 
%     % Loop over trajectory
%     for k = 1:N
%         if carriage_constraints == "centre"
%             [e_out, ~, ~] = ErrorJAtS(RailLinkage, RailCarriage, s(k), x_out(k, 1:ndof_tot)',xd_sol(k, 1:ndof_tot)');
%         elseif carriage_constraints == "both_ends"
%             [e_out, ~, ~] = ErrorJAtS_both_ends(RailLinkage, RailCarriage, s(k), x_out(k, 1:ndof_tot)',xd_sol(k, 1:ndof_tot)');
%         end
%         E_out_all(:,k) = e_out(:);
%     end
% 
%     % Time vector (assume unit timestep unless you have dt)
%     t = 1:N;
% 
%     % Plot each constraint violation
%     figure('Color','w'); hold on; grid on;
%     colors = lines(n_e);
% 
%     for i = 1:n_e
%         plot(t, E_out_all(i,:), 'LineWidth', 1.5, 'Color', colors(i,:));
%     end
% 
%     xlabel('Time step');
%     ylabel('Constraint violation (e\_out)');
%     title('Constraint Violations Over Time');
%     legend(arrayfun(@(i) sprintf('e_{out,%d}', i), 1:n_e, 'UniformOutput', false));
%     print('plots/constraint_violation_mass.png','-dpng','-r300');
% 
end

