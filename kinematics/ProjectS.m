
function [s,L] = ProjectS(Linkage,Carriage,q_b,q_c,s_prev)

% Projects carriage position onto backbone:
%
%   s = argmin ||p_backbone(s) - p_carriage||^2
%
% Inputs:
%   q_b    : backbone coordinates
%   q_c    : carriage SE(3) coordinates
%   s_prev : optional previous solution
%
% Outputs:
%   s      : projected arc-length coordinate
%   L      : total backbone length

    L = Linkage.VLinks(1).L;

    %---------------------------------------
    % Carriage position
    %---------------------------------------

    gm_end_1  = variable_expmap_g(q_c);

    %g_offset = ginv(Carriage.VLinks(1).gi);
    g_offset = Carriage.VLinks(1).gi;
    gm = gm_end_1*g_offset;  %move to centre of carriage

    p_c = gm(1:3,4);

    %---------------------------------------
    % Optional warm-start preference vs Global search
    %---------------------------------------

    if nargin >= 5 && ~isempty(s_prev)
        s = s_prev;
    else
        ds = 0.02*L;
        s_grid = linspace(0,L,20);
        %s_grid = linspace(0,L,Ngrid);
        Ngrid = length(s_grid);
        d2_best = Inf;
        s_best  = L/2;

        for k = 1:Ngrid

            T = FwdKinematicsAtS(Linkage,q_b,s_grid(k));

            p_b = T(1:3,4);

            d2 = sum((p_b - p_c).^2);

            if d2 < d2_best
                d2_best = d2;
                s_best = s_grid(k);
            end
            s = s_best;
        end
    end

    %---------------------------------------
    % Newton refinement
    %---------------------------------------

    maxIter = 30;
    tol_s   = 1e-8;

    for iter = 1:maxIter

        Tc = FwdKinematicsAtS(Linkage,q_b,s);

        p    = Tc(1:3,4);
        dpds = Tc(1:3,1);   % if x-axis is the backbone tangent

        r = p - p_c;

        g = dpds.'*r;

        if abs(g) < 1e-10
            break
        end

        %ds = -(dpds.'*r)/(dpds.'*dpds);
        denom = dpds.'*dpds + 1e-12;
        ds = -(dpds.'*r)/denom;
        
        ds = max(min(ds,0.05*L),-0.05*L); %trust region

        s_new = min(max(s + ds,0),L);

        if abs(s_new - s) < tol_s
            s = s_new;
            break
        end

        s = s_new;

    end

end
