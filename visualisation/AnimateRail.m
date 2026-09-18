function AnimateRail(Linkage, Carriage, t, x, ctrl, options)

arguments
    Linkage
    Carriage
    t
    x
    ctrl
    % Options
    options.video_name = "Dynamics"
    options.video_path = "."
    options.record = true
    % options.video_ext = ".mp4"
end

n_r = 9;

close all

% Video Settings
video_ext = ".mp4";

% Fullfile avoids conflicts with Ubuntu, since Windows use "\" and Ubuntu
% use "/".
video_file = fullfile(options.video_path, options.video_name + video_ext);

PlotParameters = Linkage.PlotParameters;



tic
tmax        = max(t);

% With fullfile, we avoid conflicts with Ubuntu/Windows
FrameRate   = PlotParameters.FrameRateValue;
if options.record
    v = VideoWriter(video_file, 'MPEG-4');
    v.FrameRate = FrameRate;
    v.Quality   = 100;
    open(v);
end

%Plot options

fh=figure(1);
%fh.Units='normalized';
fh.Units = 'pixels';
fh.Position = [100 100 1920 1080];   % Full HD
FigScale = PlotParameters.VideoResolution;
FigScale(FigScale<0.1)=0.5;
FigScale(FigScale>1)=1;
FigLocation = (1-FigScale)/2;
%fh.OuterPosition=[FigLocation FigLocation FigScale FigScale];

set(gca,'CameraPosition',PlotParameters.CameraPosition,...
    'CameraTarget',PlotParameters.CameraTarget,...
    'CameraUpVector',PlotParameters.CameraUpVector,...
    'FontSize',18)

if PlotParameters.Light
    camlight(PlotParameters.Az_light,PlotParameters.El_light)
end

axis equal
grid on
hold on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

% Set all text elements to use LaTeX interpreter
set(get(gca, 'Title'), 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'Interpreter', 'latex');
set(get(gca, 'YLabel'), 'Interpreter', 'latex');
set(get(gca, 'ZLabel'), 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

set(gca,'FontSize',12)

set(gcf, 'Renderer', 'OpenGL');

%axis limits
axis ([[-0.1 1.2*Linkage.VLinks.L] ...
        [-5*Linkage.VLinks.h{1}(0) 5*Linkage.VLinks.h{1}(0)] ...
        [-5*Linkage.VLinks.w{1}(0) 5*Linkage.VLinks.w{1}(0)]]);

drawnow


for tt=0:1/FrameRate:tmax
    cla
    %delete(findobj('type', 'patch'));
    title(strcat('t= ',num2str(tt)))
   
    xtt = interp1(t,x,tt);
    q     = xtt(1:Linkage.ndof)';
    q_mass = xtt(Linkage.ndof+1:Linkage.ndof+6)';

    %dof_start = 1;
    %g_tip    = repmat(eye(4),N,1);
    
    plotCRail(Linkage,q);
    plotLinkage(Carriage,q_mass);

    %plot carriage interface frame
    g_s = variable_expmap_g(q_mass);
    plotFrame(g_s)
    
    %plot projected tangent frame (aka where the carraiage should be)
    s_proj = ProjectS(Linkage,q,q_mass,Linkage.VLinks(1).L/2);
    p_proj = FwdKinematicsAtS(Linkage,q,s_proj);
    plotFrame(p_proj,':');

    %plot commanded carriage frame
    p_cmd = FwdKinematicsAtS(Linkage,q,ctrl.s_cmd);
    plotFrame(p_cmd, ':', 'r');

    drawnow limitrate nocallbacks;

    if options.record
        frame = getframe(gcf);
        frame.cdata = imresize(frame.cdata, [ceil(size(frame.cdata,1)/2)*2, ceil(size(frame.cdata,2)/2)*2]);
        writeVideo(v,frame);
    end
end

if options.record
    close(v);

    % Play the video only if it was recorded
    answer = questdlg('Play output video in MATLAB?','Grapical Output', ...
	    'Yes','No','Yes');
    
    if strcmp('Yes',answer)
        implay(video_file)
    end
end