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

% plot the chosen strain shape functions
%plotStrainBasis(RailLinkage)

%define properties of carriage (just a 6DOF lumped mass)

%RailCarriage = 
% RailCarriage.mass       = 105000.0; % kg
% RailCarriage.com_offset = [0.0 0.0 0.0]; % %distance from body frame
% RailCarriage.I_com     = 0.0001*eye(3); % Moment of inertia (assumes uniform & symmetrical about xyz for now)
% RailCarriage.fixed     = false;
% RailCarriage.ndof = 6;

ndof = RailLinkage.ndof+RailCarriage.ndof; %6 for the carriage

% Define Initial Conditions & Sim params
q_b0 = zeros(RailLinkage.ndof,1); % initial coordinates of rail
s0 = RailLink.L/2;%+0.2;%+1.0;

tf = 5.0; %simulation time
dt = 0.005;

tf = 0.9; %simulation time
dt = 0.005;

% calculate the initial pose of the carriage frame
g_s0 = FwdKinematicsAtS(RailLinkage,q_b0,s0); %SE(3)
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
                    s0+3,q_guess, ...
                    T_0_fixed, T_L_fixed);

[s0_new,~] = ProjectS(RailLinkage,q_equil(1:RailLinkage.ndof),q_equil(RailLinkage.ndof+1:end),s0);

plotq(RailLinkage,q_equil(1:RailLinkage.ndof))
x0 = [q_equil;zeros(ndof,1)];

%Solve DAE
[tvec_out, x_out, xdot_out] = ODEKKTSolver(RailLinkage,RailCarriage,x0,tf,dt,T_0_fixed,T_L_fixed, PD_controller );

%Plot results & animate rail

plotResults(tvec_out, x_out, xdot_out, dt, RailLinkage, RailCarriage,T_0_fixed,T_L_fixed)
AnimateRail(RailLinkage, RailCarriage, tvec_out, x_out, PD_controller)

%s_range = RailCarriage.VLinks(1).L/2:1:5.4;

s_range = 0.0:0.1:RailLinkage.VLinks(1).L;
q_equil_prev = q_guess;
q_def_z = zeros(length(s_range),1);
for i = 1:length(s_range)
    g_s = FwdKinematicsAtS(RailLinkage,q_b0,s_range(i)); %SE(3)
    q_m0 = piecewise_logmap(g_s); % carriage coordinates (in tangent space of origin)

    %define the initial system state (position and vel)
    q0 = [q_b0;q_m0];

    s = s_range(i);

    q_equil_full = FindEquilibriumFixed( ...
        RailLinkage, RailCarriage, ...
        s, q_equil_prev, ...
        T_0_fixed, T_L_fixed);

    q_carriage_se3_here = q_equil_full(end-5:end);

    T_carriage_here = variable_expmap_g(q_carriage_se3_here);

    q_def_z(i) = T_carriage_here(3,4);
    q_equil_prev = q_equil_full;
end
plot(s_range,q_def_z)
xlabel('Position along Rail (s)');
ylabel('Rail Deformation (m)');

%plotq(RailLinkage,q_equil_prev(1:RailLinkage.ndof))


[test1,test2,test3] = ErrorJAtS_both_ends(RailLinkage, RailCarriage, s0, q0, q0) ;