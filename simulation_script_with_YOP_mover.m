%Import SoRoSim Functions
run("startup.m") %run Differentiable SoRoSim startup script

%Import extra Functions
addpath('dynamics')
addpath('helper_functions')
addpath('kinematics')
addpath('visualisation')
addpath('link')
addpath('linkage')
addpath('controllers')
addpath('floating_beam')

%import properties of deformable Rail
run("RailParameters.m")
run("import_YOP_mover.m")

% plot the chosen strain shape functions
plotStrainBasis(RailLinkage)

%define properties of carriage (just a 6DOF lumped mass)
RailCarriage.mass       = 100.0; % kg
RailCarriage.com_offset = [0.0 5.0 0.0]; % %distance from body frame
RailCarriage.I_com     = 10*eye(3); % Moment of inertia (assumes uniform & symmetrical about xyz for now)
RailCarriage.fixed     = false;
RailCarriage.ndof = 6;

ndof = RailLinkage.ndof+7; %7 for the YOP payload CORRECT!!

% Define Initial Conditions & Sim params
q_b0 = zeros(RailLinkage.ndof,1); % initial coordinates of rail
s0 = RailLink.L/2;%+1.0;


tf = 2.0;%10; %simulation time
dt = 0.02;

% calculate the initial pose of the carriage frame
g_s0 = FwdKinematicsAtS(RailLinkage,q_b0,s0); %SE(3)
%q_m0_se3 = piecewise_logmap(g_s0); % carriage coordinates (in tangent space of origin)
q_m0=[g_s0(1:3,4);0;0;0;0];

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
PD_controller.kp = 1000;%0.0;%
PD_controller.kd = 1500;%0.0;%

Mover.femodel = femodel;
Mover.origins = origins;
Mover.Psi = Psi;
Mover.DOF = DOF;
Mover.m_i = m_i;
Mover.weights = weights;
Mover.rom = rom;

%Solve DAE
[tvec_out, x_out, xdot_out] = ODEKKTSolverYOPMover(RailLinkage,RailCarriage,Mover,x0,tf,dt,T_0_fixed,T_L_fixed, PD_controller );

%Plot results & animate rail
plotResultsYOP(tvec_out, x_out, xdot_out, dt, RailLinkage, RailCarriage,T_0_fixed,T_L_fixed)
AnimateRailYOP(RailLinkage, Mover, tvec_out, x_out, PD_controller,twebb,tflange)
%AnimateRailOnly(RailLinkage, tvec_out, x_out, PD_controller)

