
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
    % Global search
    %---------------------------------------

    Ngrid = 100;

    s_grid = linspace(0,L,Ngrid);

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

    end

    %---------------------------------------
    % Optional warm-start preference
    %---------------------------------------

    if nargin >= 4 && ~isempty(s_prev)

        window = 0.1*L;

        if abs(s_best - s_prev) < window
            s = s_prev;
        else
            s = s_best;
        end

    else

        s = s_best;

    end

    %---------------------------------------
    % Newton refinement
    %---------------------------------------

    maxIter = 50;
    tol_s   = 1e-8;
    tol_r   = 1e-8;

    for iter = 1:maxIter

        Tc = FwdKinematicsAtS(Linkage,q_b,s);

        p = Tc(1:3,4);

        r = p - p_c;

        if norm(r) < tol_r
            break
        end

        h = max(1e-6*L,1e-9);

        sp = min(s+h,L);
        sm = max(s-h,0);

        Tp = FwdKinematicsAtS(Linkage,q_b,sp);
        Tm = FwdKinematicsAtS(Linkage,q_b,sm);

        dpds = (Tp(1:3,4) - Tm(1:3,4)) / (sp-sm);

        denom = dpds.'*dpds;

        if denom < 1e-12
            break
        end

        ds = -(dpds.'*r)/(denom + 1e-12);

        % trust-region style clamp
        ds = max(min(ds,0.05*L),-0.05*L);

        s_new = min(max(s + ds,0),L);

        if abs(s_new - s) < tol_s
            s = s_new;
            break
        end

        s = s_new;

    end

end
