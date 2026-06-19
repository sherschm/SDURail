%Import SoRoSim Functions
run("../startup.m") %run Differentiable SoRoSim startup script

%Import extra Functions
addpath('dynamics')
addpath('helper_functions')
addpath('kinematics')
addpath('visualisation')
addpath('link')
addpath('linkage')

%import properties of deformable Rail
run("RailParameters.m")

% plot the chosen strain shape functions
plotStrainBasis(RailLinkage)

%define properties of carriage (just a 6DOF lumped mass)
RailCarriage.mass       = 10000.0; % kg
RailCarriage.com_offset = [0.0 3.0 7.0]; % %distance from body frame
RailCarriage.I_com     = 100*eye(3); % Moment of inertia (assumes uniform & symmetrical about xyz for now)
RailCarriage.fixed     = false;
RailCarriage.ndof = 6;

ndof = RailLinkage.ndof+RailCarriage.ndof; %6 for the carriage

% Define Initial Conditions & Sim params
q_b0 = zeros(RailLinkage.ndof,1); % initial coordinates of rail
s0 = RailLink.L/2+1.0;
tf = 8.0; %simulation time
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

%Solve DAE
[tvec_out, x_out, xdot_out] = ODEKKTSolver(RailLinkage,RailCarriage,x0,tf,dt,T_0_fixed,T_L_fixed);

%Plot results & animate rail
plotResults(tvec_out, x_out, dt, RailLinkage, RailCarriage,T_0_fixed,T_L_fixed)
AnimateRail(RailLinkage, tvec_out, x_out)