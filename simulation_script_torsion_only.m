run("../startup.m") %run Differentiable SoRoSim startup script
addpath('dynamics')
addpath('helper_functions')
addpath('kinematics')
addpath('visualisation')
addpath('link')
addpath('linkage')

run("RailParameters.m")

ndof = RailLinkage.ndof+6;
q_b0 = zeros(RailLinkage.ndof,1);
s0 = 2.5;
tf = 1.0;
dt = 0.001;

%Plot strain basis functions
X_vals = linspace(0,1,100);  % 100 points from 0 to 1
Phi_all = [];
for i = 1:length(X_vals)
    Phi = Phi_LegendrePolynomial(X_vals(i), RailLinkage.CVRods{1}(2).Phi_dof, RailLinkage.CVRods{1}(2).Phi_odr);
    Phi_all(:,:,i) = Phi;   % store each result
end
figure; hold on;

for j = 1:size(Phi_all,2)
    y = squeeze(Phi_all(1,j,:));
    plot(X_vals, y);
    y = squeeze(Phi_all(2,1,:));
    plot(X_vals, y);
end

xlabel('X');
ylabel('\Phi(1,j)');
title('First row of Phi vs X');
grid on;

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

%[tvec_out, x_out, xdot_out] = Baumgarte(RailLinkage,x0,tf,dt,support);

lambda0_guess = zeros(nc,1);
y0 = [x0; lambda0_guess];
ydot0 = zeros(size(y0));

[tvec_out, x_out, xdot_out, lam_out] = DAESolverTorsionOnly(RailLinkage,y0,ydot0,nc,tf,dt,support,fixed_carriage);

AnimateRail(RailLinkage, tvec_out, x_out)

plotResults(tvec_out, x_out, dt, RailLinkage)