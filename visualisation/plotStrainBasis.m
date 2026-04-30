function plotStrainBasis(RailLinkage)
    %Legendre only!
    %Plot strain basis functions
    X_vals =  linspace(0,1,100);  % 100 points from 0 to 1
    Phi_all = [];
    for i = 1:length(X_vals)
        Phi = Phi_LegendrePolynomial(X_vals(i), RailLinkage.CVRods{1}(2).Phi_dof, RailLinkage.CVRods{1}(2).Phi_odr);
        Phi_all(:,:,i) = Phi;   % store each result
    end
    figure; hold on;
    
    for j = 1:size(Phi_all,2)
        y = squeeze(Phi_all(1,j,:));
        plot(RailLinkage.VLinks.L * X_vals, y);
        y = squeeze(Phi_all(2,1,:));
        plot(RailLinkage.VLinks.L * X_vals, y);
    end
    
    xlabel('X');
    ylabel('\Phi(1,j)');
    title('First row of Phi vs X');
    grid on;

    % Save
    set(gcf, 'Units','pixels', 'Position',[100 100 1000 350]);
    set(findall(gcf,'type','axes'),'FontSize',14);
    
    set(gcf, 'PaperPositionMode','auto');
    print('plots/strain_bases.png','-dpng','-r300');
end