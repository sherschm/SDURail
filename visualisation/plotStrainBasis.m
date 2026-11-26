function plotStrainBasis(RailLinkage)

% Evaluate basis functions
X_vals = linspace(0,1,400);

Phi_all = [];
for i = 1:length(X_vals)
    Phi = RailLinkage.CVRods{1}(2).Phi_h(X_vals(i), RailLinkage.CVRods{1}(2).Phi_dof, RailLinkage.CVRods{1}(2).Phi_odr);
    Phi_torsion_only = Phi(1,:);
    Phi_torsion_only_nonzero = Phi_torsion_only(Phi_torsion_only ~= 0);
    Phi_all(:,:,i) = Phi_torsion_only_nonzero;   % store each result
end

x = RailLinkage.VLinks.L * X_vals;

% Figure
figure('Color','w',...
    'Position',[100 100 1000 300]);

ax = axes;
hold(ax,'on');

colors = lines(size(Phi_all,2));
colors(3,:) = [0.64 0.08 0.18];   % dark red

styles = { ...
    '-', ...
    '--', ...
    ':', ...
    '-.', ...
    '-', ...
    '--', ...
    ':', ...
    '-.', ...
    };
% Plot row 1 basis functions
for j = 1:size(Phi_all,2)

    y = squeeze(Phi_all(1,j,:));

    plot(x,y,...
        styles{j},'LineWidth',2,...
        'Color',colors(j,:),...
        'DisplayName',sprintf('$B_{q,%d}$',j));
end

xlabel('Position along beam, $s$ [m]', ...
    'Interpreter','latex','Color','k');

ylabel('Basis magnitude', ...
    'Interpreter','latex','Color','k');

title('Strain Basis Shapes', ...
    'Interpreter','latex','Color','k');

legend('Location','eastoutside', ...
    'Interpreter','latex');
set(legend,'Color','w')
set(legend,'TextColor','k')
set(legend,'EdgeColor','k')

grid on
box off

set(ax,...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'LineWidth',1.2,...
    'TickLabelInterpreter','latex',...
    'Color','w');
ax = gca;
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';
ax.LineWidth = 1.0;
ax.TickDir = 'out';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
grid on;
ax.GridAlpha = 0.15;
ax.MinorGridAlpha = 0.08;
ylim padded

exportgraphics(gcf,...
    'plots/strain_bases.pdf',...
    'ContentType','vector');

end