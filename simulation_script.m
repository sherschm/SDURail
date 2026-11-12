%Import SoRoSim Functions
run("startup.m") %run Differentiable SoRoSim startup script

%Import extra Functions
addpath('dynamics')
addpath('helper_functions')
addpath('kinematics')
addpath('visualisation')
addpath('link_parameters')
addpath('controllers')

%import properties of deformable Rail
run("RailParameters.m")
run("CarriageParameters.m")
%run("CarriageParameters_Offset.m")

% plot the chosen strain shape functions
%plotStrainBasis(RailLinkage)

ndof = RailLinkage.ndof+RailCarriage.ndof; %6 for the carriage

% Define Initial Conditions & Sim params
q_b0 = zeros(RailLinkage.ndof,1); % initial coordinates of rail
s0 = RailLink.L/2;%+0.2;%+1.0;

tf = 5.0; %simulation time
dt = 0.005;

tf = 1.0; %simulation time

% calculate the initial pose of the carriage frame
g_s0 = FwdKinematicsAtS(RailLinkage,q_b0,s0-RailCarriage.VLinks(1).L/2); %SE(3)
q_m0 = piecewise_logmap(g_s0); % carriage coordinates (in tangent space of origin)
   
%define the initial system state (position and vel)
q0 = [q_b0;q_m0];
v0 = zeros(ndof,1);
x0 = [q0;v0];

%Calculated rail fixed-end poses
T_full_init = FwdKinematics(RailLinkage, x0(1:RailLinkage.ndof));
T_L_fixed = T_full_init(end-3:end, :);
T_0_fixed = T_full_init(1:4,:);

%controller parameters 
PD_controller.s_cmd = 1.5; 
PD_controller.sd_cmd = 0.0;
%PD_controller.kp = 10000;
%PD_controller.kd = 15000;
PD_controller.kp = 0;
PD_controller.kd = 0;

q_guess = q0;
q_equil = FindEquilibriumFixed( ...
                    RailLinkage,RailCarriage, ...
                    s0,q_guess, ...
                    T_0_fixed, T_L_fixed);
[s0_new,~] = ProjectS(RailLinkage,RailCarriage,q_equil(1:RailLinkage.ndof),q_equil(RailLinkage.ndof+1:end),s0+RailCarriage.VLinks(1).L/2);

q0 = q_equil;
% plotq(RailLinkage,q0(1:RailLinkage.ndof))

x0 = [q0;zeros(ndof,1)];

%Solve DAE
[tvec_out, x_out, xdot_out] = ODEKKTSolver(RailLinkage,RailCarriage,x0,tf,dt,T_0_fixed,T_L_fixed, PD_controller );

%Plot results & animate rail

plotResults(tvec_out, x_out, xdot_out, dt, RailLinkage, RailCarriage,T_0_fixed,T_L_fixed)
animateSystem(RailLinkage, RailCarriage, tvec_out, x_out, PD_controller)

%s_range = RailCarriage.VLinks(1).L/2:1:5.4;
%%
s_range = RailCarriage.VLinks(1).L/2:0.2:RailLinkage.VLinks(1).L-RailCarriage.VLinks(1).L/2;%-RailCarriage.VLinks(1).L;
%s_range = 0.0:1.0:RailLinkage.VLinks(1).L;%-RailCarriage.VLinks(1).L;

q_equil_prev = q_guess;
carriage_def = zeros(length(s_range),3);
local_rotation = zeros(length(s_range),3);

for i = 1:length(s_range)
    s = s_range(i);

    q_equil_full = FindEquilibriumFixed( ...
        RailLinkage, RailCarriage, ...
        s, q_equil_prev, ...
        T_0_fixed, T_L_fixed);

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


%[test1,test2,test3] = ErrorJAtS_both_ends(RailLinkage, RailCarriage, s0, q0, q0) ;