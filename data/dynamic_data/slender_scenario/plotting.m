addpath("data\dynamic_data\slender_scenario\shape_functions_vs_nat_freq\")

range = 2:5;
frequencies = zeros(length(range),2); %first column is torsion, 2nd is bending
deformations = zeros(length(range),1); 

for i = range
    load(sprintf('%d_bases_data.mat',i))

    frequencies(i-1,1) = dominant_freq(1);
    frequencies(i-1,2) = dominant_freq(6);
    deformations(i-1) = z_carriage_equil;
end

%Bending
I = RailLink.w{1}(0)*RailLink.h{1}(0)^3/12;
EI = RailLink.E*I;
m = RailCarriage.VLinks(1).M(4,4);
f1_bending_analytical = (1/(2*pi))*sqrt(192*EI/(m*RailLink.L^3));

%
% Torsion
I_roll = RailCarriage.VLinks(1).M(1,1);
J = 0.1406* RailLink.w{1}(0)^4;
G = RailLink.G;
f1_torsion_analytical = (1/(2*pi))*sqrt(4*G*J/(I_roll*RailLink.L));

z_carriage_equil_analytical = m*9.81*RailLink.L^3/(192*EI);



%% Static deformation error plot

figure('Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% Static deformation error
ax = nexttile;

error_pct = 100*(deformations+z_carriage_equil_analytical)/z_carriage_equil_analytical;
%figure('Color','w');


ax.Color  = 'w';   % plot area
ax.XColor = 'k';   % x-axis and ticks
ax.YColor = 'k';   % y-axis and ticks
hold on 
plot(range, error_pct, ...
    '-', ...
    'Color','k', ...
    'LineWidth', 1.5);

scatter(range, error_pct, ...
    80, 'k', ...
    'filled', ...
    'MarkerEdgeColor','k');

grid on;
box on;

xticks(min(range):max(range));

xlabel('Order of Legendre basis polynomials', ...
    'Interpreter','latex', ...
    'Color','k');

ylabel('Static deflection error (\%)', ...
    'Interpreter','latex', ...
    'Color','k');
% 
% title('Convergence of static deflection', ...
%     'Interpreter','latex', ...
%     'Color','k');

set(ax,...
    'FontSize',12,...
    'LineWidth',1,...
    'GridColor',[0.8 0.8 0.8],...
    'GridAlpha',1);
ax = gca;
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.0;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;
%set(gcf,'Units','inches','Position',[1 1 5 3.5]);

ax = nexttile;

% frequency error plot
%figure('Color','w');
freq_error_pct_torsion = 100*(frequencies(:,1)-f1_torsion_analytical)/f1_torsion_analytical;
freq_error_pct_bending = 100*(frequencies(:,2)-f1_bending_analytical)/f1_bending_analytical;
for i = 1:length(range)
    freq_error_pct_torsion(i) = norm(freq_error_pct_torsion(i));
    freq_error_pct_bending(i) = norm(freq_error_pct_bending(i));
end

ax.Color  = 'w';   % plot area
ax.XColor = 'k';   % x-axis and ticks
ax.YColor = 'k';   % y-axis and ticks

hold on

plot(range, freq_error_pct_torsion, ...
    '-', ...
    'Color',[0.8500 0.3250 0.0980], ...
    'LineWidth', 1.5);

h1 = scatter(range, freq_error_pct_torsion, ...
    80, [0.8500 0.3250 0.0980], ...
    'filled', ...
    'MarkerEdgeColor', 'k');
hold on 
plot(range, freq_error_pct_bending, ...
    '-', ...
    'Color', [0 0.4470 0.7410], ...
    'LineWidth', 1.5);

h2 = scatter(range, freq_error_pct_bending, ...
    80, [0 0.4470 0.7410] , ...
    'x', ...
    'LineWidth', 1.5);

grid on;
box on;
lgd = legend([h1 h2], ...
    {'Torsional','Bending'}, ...
    'Interpreter','latex', ...
    'Location','northeast');

lgd.Box = 'off';      % optional
lgd.Color = 'w';      % white background
lgd.TextColor = 'k';  % black text
xticks(min(range):max(range));

xlabel('Order of Legendre basis polynomials', ...
    'Interpreter','latex', ...
    'Color','k');

ylabel('Natural frequency error (\%)', ...
    'Interpreter','latex', ...
    'Color','k');
% 
% title('Convergence of natural frequencies', ...
%     'Interpreter','latex', ...
%     'Color','k');

set(ax,...
    'FontSize',12,...
    'LineWidth',1,...
    'GridColor',[0.8 0.8 0.8],...
    'GridAlpha',1);
ax = gca;
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.0;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'on';
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;
set(gcf,'Units','inches','Position',[1 1 4 5.5]);
exportgraphics(gcf, ...
    'plots/shape_func_convergence.pdf', ...
    'ContentType','vector');
%%

nexttile
ax = axes;
hold(ax,'on')

% White axes background
set(ax,...
    'Color','w',...
    'XColor','k',...
    'YColor','k',...
    'FontName','Times',...
    'FontSize',11,...
    'LineWidth',1.0);

grid on
box on

% Responses
plot(f, amplitude(:,1)./max(abs(amplitude(:,1))), ...
    'LineWidth',1.3,'Color', [0.8500 0.3250 0.0980]);

plot(f, amplitude(:,6)./max(abs(amplitude(:,6))), ...
    'LineWidth',1.3,'Color', [0 0.4470 0.7410]);

% Analytical frequencies
xline(f1_torsion_analytical,'--','LineWidth',1.8,'Color', [0.8500 0.3250 0.0980]);
xline(f1_bending_analytical,'--','LineWidth',1.8,'Color', [0 0.4470 0.7410]);

xlabel('Frequency (Hz)', ...
    'FontSize',12, ...
    'Color','k');

ylabel('Normalised Amplitude (-)', ...
    'FontSize',12, ...
    'Color','k');

legend({'Roll response',...
    'Vertical response',...
    'Analytical torsion',...
    'Analytical bending'}, ...
    'Location','best');

set(legend,'Color','w')
set(legend,'TextColor','k')
set(legend,'EdgeColor','k')

% Squish vertically
ylim([0 1.1])      % choose range of interest
xlim([0 150])      % choose range of interest

% ylim([0 0.25])   % even flatter
% ylim([0 0.15])   % highlights small peaks

% Make figure itself shorter
set(gcf,'Position',[100 100 700 250]);

set(gcf,'Renderer','painters');

exportgraphics(gcf, ...
    'plots/carriage_FFT.pdf', ...
    'ContentType','vector');

%Error
