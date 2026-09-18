%Import SoRoSim Functions
run("startup.m") %run Differentiable SoRoSim startup script

%Import extra Functions
addpath('dynamics')
addpath('helper_functions')
addpath('kinematics')
addpath('visualisation')
addpath('link')
addpath('controllers')

%import properties of deformable Rail
run("RailParameters.m")
run("CarriageParameters.m")

% plot the chosen strain shape functions
plotStrainBasis(RailLinkage)

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
s0 = RailLink.L/2+0.2;%+1.0;

tf = 5.0; %simulation time
dt = 0.005;

tf = 1; %simulation time
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
q_equil = FindEquilibrium( ...
                    RailLinkage,RailCarriage, ...
                    s0,q_guess, ...
                    T_0_fixed, T_L_fixed)

plotq(RailLinkage,q_equil(1:RailLinkage.ndof))
x0 = [q_equil;zeros(ndof,1)];

%Solve DAE
[tvec_out, x_out, xdot_out] = ODEKKTSolver(RailLinkage,RailCarriage,x0,tf,dt,T_0_fixed,T_L_fixed, PD_controller );

%Plot results & animate rail

plotResults(tvec_out, x_out, xdot_out, dt, RailLinkage, RailCarriage,T_0_fixed,T_L_fixed)
AnimateRail(RailLinkage, RailCarriage, tvec_out, x_out, PD_controller)
