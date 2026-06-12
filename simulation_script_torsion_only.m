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

%define properties of carriage
RailCarriage.mass       = 0.01;%0000.0; % kg
RailCarriage.com_offset = [0.0 0.0 0.0]; % %distance from body frame
RailCarriage.I     = 1*eye(3); % Moment of inertia (assumes uniform & symmetrical about xyz for now)
RailCarriage.fixed     = false;
RailCarriage.ndof = 6;

ndof = RailLinkage.ndof+RailCarriage.ndof; %6 for the carriage
q_b0 = zeros(RailLinkage.ndof,1); % initial coordinates of rail
%s0 = 2.5; % initial position of carriage along rail
s0 = RailLink.L/2+1.0;
tf = 10.0; %simulation time
dt = 0.001;

% calculate the initial pose of the carriage frame
g_s0 = FwdKinematicsAtS(RailLinkage,q_b0,s0); %SE(3)
q_m0 = piecewise_logmap(g_s0); % carriage coordinates (in tangent space of origin)
   
%Jr = SE3RightJacobianFromPose(g);
%define the initial system state (position and vel)
q0 = [q_b0;q_m0];
v0 = zeros(ndof,1);
x0 = [q0;v0];

%Solve DAE
%[tvec_out, x_out, xdot_out, lam_out] = DAESolverTorsionOnly(RailLinkage,RailCarriage,x0,tf,dt);
%[tvec_out, x_out, xdot_out, lam_out] =DAESolverImplicit(RailLinkage,x0,tf,dt,RailCarriage);

[tvec_out, x_out, xdot_out] = ODESolverTorsionOnly(RailLinkage,RailCarriage,x0,tf,dt);

%Plot results & animate rail
plotResults(tvec_out, x_out, dt, RailLinkage)
AnimateRail(RailLinkage, tvec_out, x_out)

 %plotq(RailLinkage,x_out(3000,1:RailLinkage.ndof))
 %plotFrame(FwdKinematicsAtS(RailLinkage,x_out(3000,1:RailLinkage.ndof),2.7))
% plotRail(RailLinkage,x_out(end,1:RailLinkage.ndof))
% FwdKinematicsAtS(RailLinkage,x_out(end,:),s0)

%[ID_no_Mqdd, tau, e, dID_dq, Cb, Mb] = DAEJacobians(RailLinkage,0.0,q_b0,zeros(RailLinkage.ndof,1),zeros(RailLinkage.ndof,1),[],[],[]);
%GeneralizedMassMatrix_corrected(RailLinkage,q_b0)