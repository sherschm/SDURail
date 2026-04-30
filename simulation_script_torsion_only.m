run("../startup.m") %run Differentiable SoRoSim startup script
addpath('dynamics')
addpath('helper_functions')
addpath('kinematics')
addpath('visualisation')
addpath('link')
addpath('linkage')

run("RailParameters.m")

plotStrainBasis(RailLinkage)

ndof = RailLinkage.ndof+6;
q_b0 = zeros(RailLinkage.ndof,1);
s0 = 2.5;
tf = 5;%1.2; %
dt = 0.001;

g_s0 = FwdKinematicsAtS(RailLinkage,q_b0,s0);
q_m0 = piecewise_logmap(g_s0);

x0 = [q_b0;q_m0;zeros(ndof,1)];

fixed_carriage = false;

if fixed_carriage == true
    nc_carriage = 6;
else
    nc_carriage = 5;
end

%support = "pin-pin"; 
%nc = nc_carriage+ 6;

%support = "pin-fixed"; 
%nc = nc_carriage+ 9;

support = "fixed-fixed"; 
%nc = nc_carriage + 12;
nc = nc_carriage + 7;
%nc = nc_carriage + 6;

%[tvec_out, x_out, xdot_out] = Baumgarte(RailLinkage,x0,tf,dt,support);

lambda0_guess = zeros(nc,1);
y0 = [x0; lambda0_guess];
ydot0 = zeros(size(y0));

RailCarriage.mass       = 80000.0; % kg
RailCarriage.com_offset = [0 1.5 0]; % %distance from body frame
RailCarriage.I     = 0.0001; % Moment of inertia (assumes uniform & symmetrical about xyz for now)
RailCarriage.fixed     = false;


[tvec_out, x_out, xdot_out, lam_out] = DAESolverTorsionOnly(RailLinkage,y0,ydot0,nc,tf,dt,support,RailCarriage);

plotResults(tvec_out, x_out, dt, RailLinkage)

AnimateRail(RailLinkage, tvec_out, x_out)

% plotq(RailLinkage,x_out(end,1:RailLinkage.ndof))
% plotRail(RailLinkage,x_out(end,1:RailLinkage.ndof))
% 
% FwdKinematicsAtS(RailLinkage,x_out(end,:),s0)