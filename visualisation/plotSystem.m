function plotSystem(Linkage, Carriage, q, q_mass, ctrl)

close all

PlotParameters = Linkage.PlotParameters;

figure(1);
clf

set(gca,...
    'CameraPosition', PlotParameters.CameraPosition,...
    'CameraTarget', PlotParameters.CameraTarget,...
    'CameraUpVector', PlotParameters.CameraUpVector,...
    'FontSize', 18)

if PlotParameters.Light
    camlight(PlotParameters.Az_light, ...
        PlotParameters.El_light)
end

axis equal
grid on
hold on

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

set(get(gca,'Title'),'Interpreter','latex');
set(get(gca,'XLabel'),'Interpreter','latex');
set(get(gca,'YLabel'),'Interpreter','latex');
set(get(gca,'ZLabel'),'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex');

set(gcf,'Renderer','OpenGL');

axis([ ...
    -0.1 1.2*Linkage.VLinks.L ...
    -2*Linkage.VLinks.h{1}(0) 2*Linkage.VLinks.h{1}(0) ...
    -2*Linkage.VLinks.w{1}(0) 2*Linkage.VLinks.w{1}(0)]);

title('Rail configuration')

% Rail
plotCRail(Linkage,q);

% Carriage
plotLinkage(Carriage,q_mass);

% Actual carriage frame
g_s_end = variable_expmap_g(q_mass);
gi = Carriage.VLinks(1).gi;
g_s= g_s_end*gi;
plotFrame(g_s)


% Projected frame on rail
s_proj = ProjectS( ...
    Linkage,...
    Carriage,...
    q,...
    q_mass,...
    Linkage.VLinks(1).L/2);

p_proj = FwdKinematicsAtS(Linkage,q,s_proj);
plotFrame(p_proj,':');

% Commanded frame
if nargin > 4 && ~isempty(ctrl)
    p_cmd = FwdKinematicsAtS(Linkage,q,ctrl.s_cmd);
    plotFrame(p_cmd,':','r');
end

drawnow

end