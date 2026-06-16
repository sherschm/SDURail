% function [s,L] = ProjectS(Linkage, q_b, q_c,s_prev)
% 
%     % ----------------------------
%     % parameters
%     % ----------------------------
%     tol     = 1e-9;
%     maxiter = 5;
%     reg     = 1e-10;
% 
%     L = Linkage.VLinks(1).L;
% 
%     % initial guess
%     s = L/2;
% 
%     % ----------------------------
%     % carriage position
%     % ----------------------------
%     gm  = variable_expmap_g(q_c);   % 4x4 pose
%     p_c = gm(1:3,4);
% 
%     % ----------------------------
%     % helper: beam position at s
%     % ----------------------------
%     function pb = get_beam_pos(Linkage, q_b, sigma)
%         T = FwdKinematicsAtS(Linkage, q_b, sigma);
%         pb = T(1:3,4);
%     end
% 
%     fk_pos = @(sigma) get_beam_pos(Linkage,q_b,sigma);
% 
%     % ----------------------------
%     % scalar objective
%     % Phi(s) = 0.5 * ||pb - pc||^2
%     % ----------------------------
%     Phi = @(sigma) 0.5 * sum((fk_pos(sigma) - p_c).^2);
% 
%     % ----------------------------
%     % Newton iteration
%     % ----------------------------
%     for iter = 1:maxiter
% 
%         % first derivative
%         d1 = derivative_fd(Phi, s);
% 
%         % convergence check
%         if abs(d1) < tol
%             return;
%         end
% 
%         % second derivative
%         d2 = derivative_fd(@(sigma) derivative_fd(Phi, sigma), s);
% 
%         % regularisation
%         if abs(d2) < reg
%             d2 = d2 + sign(d2 + eps) * reg;
%         end
% 
%         % Newton step
%         ds = -d1 / d2;
%         s_new = s + ds;
% 
%         % optional clamp (recommended for stability)
%         % s_new = max(0, min(L, s_new));
% 
%         % stopping condition
%         if abs(s_new - s) < 1e-14
%             s = s_new;
%             fprintf('proj iter = %d\n', iter);
%             break;
%         end
% 
%         s = s_new;
%     end
% 
% end
% 
function [s,L] = ProjectS(Linkage, q_b, q_c, s_prev)
    maxIter = 5;

    L = Linkage.VLinks(1).L;

    % warm start
    if nargin < 4 || isempty(s_prev)
        s = L/2;
    else
        s = s_prev;
    end

    % carriage position
    gm  = variable_expmap_g(q_c);
    p_c = gm(1:3,4);

    for k = 1:maxIter

        T = FwdKinematicsAtS(Linkage, q_b, s);
        p_b = T(1:3,4);
        t = T(1:3,1);

        r = p_b - p_c;
        ds = - dot(r,t) / dot(t,t);

        s = s + 0.7*ds;
        s = max(0,min(L,s));

       % if abs(g) < tol
%             converged = true;
%             break;
%         end
    end

end

% function [s,info] = ProjectS(Linkage,q_b,q_c,s0)
% 
%     L = Linkage.VLinks(1).L;
% 
%     tol = 1e-10;
%     maxIter = 1000;
% 
%     s = min(max(s0,0),L);
% 
%     gm  = variable_expmap_g(q_c);
%     p_c = gm(1:3,4);
% 
%     converged = false;
% 
%     for k = 1:maxIter
% 
%         %[p_b,t,ts] = BeamPointAndDerivatives(Linkage,q_b,s);
%         T_r = FwdKinematicsAtS(Linkage, q_b, s);
%         p_r = T_r(1:3,4);
%         %J_r = JacobianAtS(Linkage, q_b, s);
%         t_full = dinamico_Adjoint(T_r)*dfzcds(Linkage,q_b,s);
%         t = t_full(4:6);
% 
%         r = p_r - p_c;
% 
%         g  = dot(r,t);
% 
%         gs = dot(t,t); %+ dot(r,ts); %drop curvature term (add later?)
% 
%         if abs(g) < tol
%             converged = true;
%             break;
%         end
% 
%         ds = -g/gs;
% 
%         s = s + ds;
% 
%         s = min(max(s,0),L);
%     end
% 
%     info.g = g;
%     info.gs = gs;
%     info.iterations = k;
%     info.converged = converged;
%     disp(k)
%     disp(converged)
% end