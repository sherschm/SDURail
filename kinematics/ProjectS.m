% function s = ProjectS(Linkage, q_b, q_c)
% 
% % ----------------------------
% % parameters
% % ----------------------------
% tol     = 1e-9;
% maxiter = 100;
% reg     = 1e-10;
% 
% L = Linkage.VLinks(1).L;
% 
% % initial guess
% s = L/2;
% 
% % ----------------------------
% % carriage position
% % ----------------------------
% gm  = variable_expmap_g(q_c);   % 4x4 pose
% p_c = gm(1:3,4);
% 
% % ----------------------------
% % helper: beam position at s
% % ----------------------------
% function pb = get_beam_pos(Linkage, q_b, sigma)
%     T = FwdKinematicsAtS(Linkage, q_b, sigma);
%     pb = T(1:3,4);
% end
% 
% fk_pos = @(sigma) get_beam_pos(Linkage,q_b,sigma);
% 
% % ----------------------------
% % scalar objective
% % Phi(s) = 0.5 * ||pb - pc||^2
% % ----------------------------
% Phi = @(sigma) 0.5 * sum((fk_pos(sigma) - p_c).^2);
% 
% % ----------------------------
% % Newton iteration
% % ----------------------------
% for iter = 1:maxiter
% 
%     % first derivative
%     d1 = derivative_fd(Phi, s);
% 
%     % convergence check
%     if abs(d1) < tol
%         return;
%     end
% 
%     % second derivative
%     d2 = derivative_fd(@(sigma) derivative_fd(Phi, sigma), s);
% 
%     % regularisation
%     if abs(d2) < reg
%         d2 = d2 + sign(d2 + eps) * reg;
%     end
% 
%     % Newton step
%     ds = -d1 / d2;
%     s_new = s + ds;
% 
%     % optional clamp (recommended for stability)
%     % s_new = max(0, min(L, s_new));
% 
%     % stopping condition
%     if abs(s_new - s) < 1e-14
%         s = s_new;
%         fprintf('proj iter = %d\n', iter);
%         break;
%     end
% 
%     s = s_new;
% end
% 
% end

function s = ProjectS(Linkage, q_b, q_c, s_prev)

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
    
    % % beam position
    % T = FwdKinematicsAtS(Linkage, q_b, s);
    % p_b = T(1:3,4);
    % 
    % % tangent direction (WORLD FRAME)
    % t = T(1:3,1);   % <-- THIS is the key fix
    % 
    % % projection step
    % r = p_b - p_c;
    % 
    % denom = dot(t,t);
    % if denom < 1e-12
    %     return;
    % end
    % 
    % ds = - dot(r,t) / denom;
    % 
    % % damping
    % alpha = 0.7;
    % s = s + alpha * ds;
    % 
    % % clamp
    % s = max(0, min(L, s));
    
    for k = 1:2
        T = FwdKinematicsAtS(Linkage, q_b, s);
        p_b = T(1:3,4);
        t = T(1:3,1);
    
        r = p_b - p_c;
        ds = - dot(r,t) / dot(t,t);
    
        s = s + 0.7*ds;
        s = max(0,min(L,s));
    end

end