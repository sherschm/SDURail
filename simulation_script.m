%Import SoRoSim Functions
run("startup.m") %run Differentiable SoRoSim startup script

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
RailCarriage.com_offset = [0.0 5.0 0.0]; % %distance from body frame
RailCarriage.I_com     = 10*eye(3); % Moment of inertia (assumes uniform & symmetrical about xyz for now)
RailCarriage.fixed     = false;
RailCarriage.ndof = 6;

ndof = RailLinkage.ndof+RailCarriage.ndof; %6 for the carriage

% Define Initial Conditions & Sim params
q_b0 = zeros(RailLinkage.ndof,1); % initial coordinates of rail
s0 = RailLink.L/2;%+1.0;
tf = 100; %simulation time
dt = 0.02;
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
plotResults(tvec_out, x_out, xdot_out, dt, RailLinkage, RailCarriage,T_0_fixed,T_L_fixed)
AnimateRail(RailLinkage, tvec_out, x_out)


%AnimateRail(RailLinkage, tvec_out(end-1:end), x_out(end-1:end,:))

%test 
% 
q = 0.5*ones(RailLinkage.ndof,1);
qd = 0.5*ones(RailLinkage.ndof,1);
qdd = 0.5*ones(RailLinkage.ndof,1);
Jd_full = Jacobiandot(RailLinkage, q, qd);
Jd_full_test = JacobiandotAtS(RailLinkage, q, qd,RailLinkage.VLinks.L);

% [ID_no_Mqdd, tau, e, dID_dq, Cb, Mb] = DAEJacobians(RailLinkage,0.0,q,qd,zeros(RailLinkage.ndof,1),[],[],[]);
% 
% [ID_no_MqddCqd, tau, e, dID_dq, Cb, Mb] = DAEJacobians(RailLinkage,0.0,q,zeros(RailLinkage.ndof,1),zeros(RailLinkage.ndof,1),[],[],[]);
% 
% [ID, tau, e, dID_dq, Cb, Mb] = DAEJacobians(RailLinkage,0.0,q,qd,qdd,[],[],[]);
% 
% ID_test = ID_no_Mqdd + Mb*qdd;
% 
% ID