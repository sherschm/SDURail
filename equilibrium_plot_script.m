%%
%s_range = 0.0:1.0:RailLinkage.VLinks(1).L;%-RailCarriage.VLinks(1).L;

%%Import benchmark data.
addpath(".\data\equil_data\")
plot_script2_config1

constraint_configs = ["centre"; "both_ends"];
%constraint_configs = ["centre"; "both_ends";"three_points"];
names = ["Constraint (a)"; "Constraint (b)"; "Constraint (c)"];

s_here = RailLink.L/2;%+RailCarriage.VLinks(1).L/2;
%s_here = RailLink.L/2;

s_range = 0.0:0.001:RailLinkage.VLinks(1).L;%-RailCarriage.VLinks(1).L;

local_rotation = zeros(length(s_range),3,length(constraint_configs));
carriage_def   = zeros(length(s_range),3,length(constraint_configs));

for j = 1:length(constraint_configs)
    q_equil_full = FindEquilibriumFixed(RailLinkage, RailCarriage, ...
        s_here, q_equil, ...
        T_0_fixed, T_L_fixed,constraint_configs(j));

    for i = 1:length(s_range)
        s = s_range(i);
    
        q_carriage_se3_here = q_equil_full(end-5:end);
    
        T_carriage_here_end = variable_expmap_g(q_carriage_se3_here);
        
        T_carriage_here= FwdKinematicsAtS(RailLinkage,q_equil_full(1:RailLinkage.ndof),s);
    
        local_rotation(i,:,j) = LogSO3(T_carriage_here(1:3,1:3));
    
        carriage_def(i,:,j) = T_carriage_here(1:3,4);
    
    end
end

%% Lateral plots
figure('Color','w');

c_y = [0.0000 0.4470 0.7410];   % blue
c_green = [0.0000 0.6000 0.3000];

t = tiledlayout(1,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% Vertical deformation
nexttile;
hold on 
cols = lines(length(constraint_configs));

styles = {'-','--',':','-.'};   % add more if needed

h = [];
leg = {};

for j = 1:length(constraint_configs)

    % hy = plot(s_range, carriage_def(:,2,j), ...
    %     styles{j}, ...
    %     'LineWidth',2.0, ...
    %     'Color',c_y);
    hold on

    hz = plot(s_range, carriage_def(:,3,j), ...
        styles{j}, ...
        'LineWidth',2.0);
    hold on
    %h = [h hy hz];
    h = [h hz];
    %leg{end+1} = sprintf('$y$, %s', names{j});
    leg{end+1} = sprintf(names{j});
end
h_ansys = plot(z_data(:,1), z_data(:,2), ...
    'o', ...
    'MarkerSize', 4, ...
    'Color', c_green);

h = [h h_ansys];
leg{end+1} = sprintf('ANSYS', names{j});

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('Vertical rail deformation (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);


legend(h, leg, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex', TextColor ='black');
title('(a)','Interpreter','latex','Interpreter','latex', ...
    'Color','k');
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


% Angular deformation
h = [];
nexttile;


for j = 1:length(constraint_configs)

    % h_pitch = plot(s_range, -local_rotation(:,2,j), ...
    %     styles{j}, ...
    %     'LineWidth',2.0);


    % Smooth the profile
    z_smooth = smoothdata(carriage_def(:,3,j),'sgolay',501);
    theta = gradient(z_smooth,s_range);

    local_rotation_numerical = atan(theta); % radians
    h_pitch = plot(s_range, local_rotation_numerical, ...
            styles{j}, ...
         'LineWidth',2.0);

    h = [h h_pitch];

    hold on
    % plot(s_range, -local_rotation(:,3,j), ...
    %     ':', ...
    %     'LineWidth',2.5, ...
    %     'Color',[0.9290 0.6940 0.1250]);

    %leg{end+1} = sprintf('$y$, %s', names{j});
    leg{end+1} = sprintf(names{j});
end
% ansys_local_rotation = zeros(size(z_data,1)-1,1);
% for i = 1:size(z_data,1)-1
%     ansys_local_rotation(i) = (z_data(i+1,2)-z_data(i,2))/(z_data(i+1,1)-z_data(i,1));
% end

x = z_data(:,1);
z = z_data(:,2);

% Smooth the profile
%z_smooth = smoothdata(z,'sgolay',21);
z_smooth = smoothdata(z,'sgolay',21);
theta = gradient(z_smooth,x);

ansys_local_rotation = atan(theta); % radians

h_ansys = plot(z_data(:,1), ansys_local_rotation, ...
    'o', ...
    'MarkerSize', 4, ...
    'Color', c_green);

h = [h h_ansys];

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('Angular Pitch Deformation (rad)', ...
    'Interpreter','latex', ...
    'FontSize',15);

legend(h, leg, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex', TextColor ='black');
title('(b)','Interpreter','latex','Interpreter','latex', ...
    'Color','k');
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

% Figure size and export
set(gcf,'Units','inches','Position',[1 1 10 3.5]);

exportgraphics(gcf, ...
    'rail_carriage_deformation_planar.pdf', ...
    'ContentType','vector');

%plotSystem(RailLinkage,RailCarriage,q_equil_full(1:RailLinkage.ndof),q_equil_full(RailLinkage.ndof+1:end))

%% Torsional PLots:

config2;
figure('Color','w');

c_y = [0.0000 0.4470 0.7410];   % blue
c_green = [0.0000 0.6000 0.3000];

t = tiledlayout(3,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

nexttile;
hold on 
cols = lines(length(constraint_configs));

styles = {'-','--',':','-.'};   % add more if needed

% x-def
h = [];
leg = {};

for j = 1:length(constraint_configs)

    % hy = plot(s_range, carriage_def(:,2,j), ...
    %     styles{j}, ...
    %     'LineWidth',2.0, ...
    %     'Color',c_y);
    hold on

    hy = plot(s_range, carriage_def(:,2,j), ...
        styles{j}, ...
        'LineWidth',2.0);
    hold on
    %h = [h hy hz];
    h = [h hy];
    %leg{end+1} = sprintf('$y$, %s', names{j});
    leg{end+1} = sprintf(names{j});
end
% h_ansys = plot(x_data(:,1), -x_data(:,2), ...
%     'o', ...
%     'MarkerSize', 4, ...
%     'Color', c_green);
% 
% h = [h h_ansys];
% leg{end+1} = sprintf('ANSYS', names{j});

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('Lateral $y$ rail deformation (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);


legend(h, leg, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex', TextColor ='black');
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

% z-def
nexttile;
hold on 
h = [];
leg = {};

for j = 1:length(constraint_configs)

    % hy = plot(s_range, carriage_def(:,2,j), ...
    %     styles{j}, ...
    %     'LineWidth',2.0, ...
    %     'Color',c_y);
    hold on

    hz = plot(s_range, carriage_def(:,3,j), ...
        styles{j}, ...
        'LineWidth',2.0);
    hold on
    %h = [h hy hz];
    h = [h hz];
    %leg{end+1} = sprintf('$y$, %s', names{j});
    leg{end+1} = sprintf(names{j});
end

mean_z_def_ansys = (z_r_data(:,2)+z_l_data(:,2))/2;

h_ansys = plot(z_r_data(:,1), mean_z_def_ansys, ...
    'o', ...
    'MarkerSize', 4, ...
    'Color', c_green);
h = [h h_ansys];
leg{end+1} = 'ANSYS';

% h_ansys = plot(z_data(:,1), z_data(:,2), ...
%     'o', ...
%     'MarkerSize', 4, ...
%     'Color', c_green);
% 
% h = [h h_ansys];
%leg{end+1} = sprintf('ANSYS', names{j});

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('Vertical rail deformation (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);


legend(h, leg, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex', TextColor ='black');
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

% Angular deformation
h = [];
nexttile;
for j = 1:length(constraint_configs)
    % plot(s_range, local_rotation(:,1,j), ...
    %     '-', ...
    %     'LineWidth',2.0, ...
    %     'Color',[0 0.4470 0.7410]);
    % hold on
    
    h_pitch = plot(s_range, -local_rotation(:,3,j), ...
        styles{j}, ...
        'LineWidth',2.0);
    h = [h h_pitch];

    hold on
    % plot(s_range, -local_rotation(:,3,j), ...
    %     ':', ...
    %     'LineWidth',2.5, ...
    %     'Color',[0.9290 0.6940 0.1250]);

    %leg{end+1} = sprintf('$y$, %s', names{j});
    leg{end+1} = sprintf(names{j});
end
% ansys_local_rotation = zeros(size(z_data,1)-1,1);
% for i = 1:size(z_data,1)-1
%     ansys_local_rotation(i) = (z_data(i+1,2)-z_data(i,2))/(z_data(i+1,1)-z_data(i,1));
% end
% 
% x = z_data(:,1);
% z = z_data(:,2);
% 
% % Smooth the profile
% z_smooth = smoothdata(z,'sgolay',21);
% theta = gradient(z_smooth,x);
% 
% ansys_local_rotation = atan(theta); % radians
% 
% h_ansys = plot(z_data(:,1), ansys_local_rotation, ...
%     'o', ...
%     'MarkerSize', 4, ...
%     'Color', c_green);
% 
% h = [h h_ansys];

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('Angular Yaw Deformation (rad)', ...
    'Interpreter','latex', ...
    'FontSize',15);

legend(h, leg, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex', TextColor ='black');
title('(b)','Interpreter','latex','Interpreter','latex', ...
    'Color','k');
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

% Roll deformation
h = [];
nexttile;
for j = 1:length(constraint_configs)
    % plot(s_range, local_rotation(:,1,j), ...
    %     '-', ...
    %     'LineWidth',2.0, ...
    %     'Color',[0 0.4470 0.7410]);
    % hold on

    h_roll = plot(s_range, local_rotation(:,1,j), ...
        styles{j}, ...
        'LineWidth',2.0);
    h = [h h_roll];

    hold on
    % plot(s_range, -local_rotation(:,3,j), ...
    %     ':', ...
    %     'LineWidth',2.5, ...
    %     'Color',[0.9290 0.6940 0.1250]);

    %leg{end+1} = sprintf('$y$, %s', names{j});
    leg{end+1} = sprintf(names{j});
end

delta_ansys = z_r_data(:,2)- z_l_data(:,2);
rail_width = RailLinkage.VLinks(1).h{1}(0);
roll_ansys = zeros(length(delta_ansys),1);
for i = 1:length(delta_ansys)
    roll_ansys(i) = -asin(delta_ansys(i)/rail_width);
end

h_ansys = plot(z_r_data(:,1), roll_ansys, ...
    'o', ...
    'MarkerSize', 4, ...
    'Color', c_green);
leg{end+1} = 'ANSYS';

h = [h h_ansys];

xlabel('Position Along Rail, $s$ (m)', ...
    'Interpreter','latex', ...
    'FontSize',15);

ylabel('Angular Roll Deformation (rad)', ...
    'Interpreter','latex', ...
    'FontSize',15);

legend(h, leg, ...
    'Location','best', ...
    'Box','off', ...
    'Interpreter','latex', TextColor ='black');
title('(b)','Interpreter','latex','Interpreter','latex', ...
    'Color','k');
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
% Figure size and export
set(gcf,'Units','inches','Position',[1 1 10 7]);

exportgraphics(gcf, ...
    'rail_carriage_deformation_torsion.pdf', ...
    'ContentType','vector');

%plotSystem(RailLinkage,RailCarriage,q_equil_full(1:RailLinkage.ndof),q_equil_full(RailLinkage.ndof+1:end))



%% sliding carriage centrepoint defomration plot


s_range = RailCarriage.VLinks(1).L/2:0.24:RailLinkage.VLinks(1).L-RailCarriage.VLinks(1).L/2;%-RailCarriage.VLinks(1).L;

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