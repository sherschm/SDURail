%%
s_range = RailCarriage.VLinks(1).L/2:0.24:RailLinkage.VLinks(1).L-RailCarriage.VLinks(1).L/2;%-RailCarriage.VLinks(1).L;
%s_range = 0.0:1.0:RailLinkage.VLinks(1).L;%-RailCarriage.VLinks(1).L;

q_equil_prev = q_guess;
carriage_def = zeros(length(s_range),3);
local_rotation = zeros(length(s_range),3);

for i = 1:length(s_range)
    s = s_range(i);

    q_equil_full = FindEquilibriumFixed( ...
        RailLinkage, RailCarriage, ...
        s, q_equil_prev, ...
        T_0_fixed, T_L_fixed,carriage_constraints);

    q_carriage_se3_here = q_equil_full(end-5:end);

    T_carriage_here_end = variable_expmap_g(q_carriage_se3_here);
    
     %gi = RailCarriage.VLinks(1).gi;
     %T_carriage_here= T_carriage_here_end*gi;

     T_carriage_here= FwdKinematicsAtS(RailLinkage,q_equil_full(1:RailLinkage.ndof),s);
    %T_carriage_here= T_carriage_here_end;
    local_rotation(i,:) = LogSO3(T_carriage_here(1:3,1:3));

    carriage_def(i,:) = T_carriage_here(1:3,4);
    q_equil_prev = q_equil_full;
    plotSystem(RailLinkage,RailCarriage,q_equil_prev(1:RailLinkage.ndof),q_equil_prev(RailLinkage.ndof+1:end))

end

%s_range = s_range+RailCarriage.VLinks(1).L/2;%+RailCarriage.VLinks(1).L/2;

%Move deformation computation to the centre of the carriage:


figure('Color','w');

t = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Vertical deformation
nexttile;

plot(s_range, carriage_def(:,2), ...
    '--', ...
    'LineWidth',2.0, ...
    'Color',[0.8500 0.3250 0.0980]);
hold on

plot(s_range, carriage_def(:,3), ...
    ':', ...
    'LineWidth',2.5, ...
    'Color',[0.9290 0.6940 0.1250]);

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',12);

ylabel('Deformation (m)', ...
    'Interpreter','latex', ...
    'FontSize',12);

legend({'$y$','$z$'}, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex');
title('(a)','Interpreter','latex')
ax = gca;
ax.FontSize = 11;
ax.LineWidth = 1.0;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;

% Angular deformation
nexttile;

plot(s_range, local_rotation(:,1), ...
    '-', ...
    'LineWidth',2.0, ...
    'Color',[0 0.4470 0.7410]);
hold on

plot(s_range, local_rotation(:,2), ...
    '--', ...
    'LineWidth',2.0, ...
    'Color',[0.8500 0.3250 0.0980]);

plot(s_range, local_rotation(:,3), ...
    ':', ...
    'LineWidth',2.5, ...
    'Color',[0.9290 0.6940 0.1250]);

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',12);

ylabel('Angular Deformation (rad)', ...
    'Interpreter','latex', ...
    'FontSize',12);

legend({'Roll','Pitch','Yaw'}, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex');
title('(b)','Interpreter','latex')
ax = gca;
ax.FontSize = 11;
ax.LineWidth = 1.0;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;

% Figure size and export
set(gcf,'Units','inches','Position',[1 1 10 2.5]);

exportgraphics(gcf, ...
    'rail_carriage_deformation.pdf', ...
    'ContentType','vector');