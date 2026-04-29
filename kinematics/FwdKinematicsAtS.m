%Function that calculates the forward kinematics of the linkage from the base to the linkage tip at every significant points
%Last modified by Anup Teejo Mathew 29.11.2024
function g = FwdKinematicsAtS(Linkage,q,s) %i_here is link index, j_here is division (0 for joint)
% Forward kinematics evaluated up to arc-length position s ∈ [0,1]
% If s is omitted → full forward kinematics is returned

if isrow(q)
    q = q';
end

N = Linkage.N;
g_ini = Linkage.g_ini;
g_tip = repmat(eye(4), N, 1);
iLpre = Linkage.iLpre;

g_here = eye(4);

dof_start = 1;
s_acc = 0;

for i = 1:N

    % -------------------------
    % Base transform / previous link
    % -------------------------
    if iLpre(i) > 0
        g_here = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:) * ...
                 g_ini((i-1)*4+1:i*4,:);
    else
        g_here = g_ini((i-1)*4+1:i*4,:);
    end

    % -------------------------
    % Joint
    % -------------------------
    dof_here = Linkage.CVRods{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    Phi_here = Linkage.CVRods{i}(1).Phi;
    xi_star  = Linkage.CVRods{i}(1).xi_star;

    if dof_here == 0
        g_joint = eye(4);
    else
        xi = Phi_here*q_here + xi_star;
        g_joint = variable_expmap_g(xi);
    end

    g_here = g_here * g_joint;

    % check early stop
    if s <= s_acc
        g = g_here;
        return
    end

    % -------------------------
    % Rigid link (if any)
    % -------------------------
    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype == 'r'

        gi = Linkage.VLinks(Linkage.LinkIndex(i)).gi;
        g_here = g_here * gi;

        if ~full_eval && s <= s_acc
            g = g_here;
            return
        end

        gf = Linkage.VLinks(Linkage.LinkIndex(i)).gf;
        g_here = g_here * gf;
    end

    dof_start = dof_start + dof_here;

    % -------------------------
    % Flexible pieces
    % -------------------------
    for j = 1:Linkage.VLinks(Linkage.LinkIndex(i)).npie - 1

        dof_here = Linkage.CVRods{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);

        xi_star = Linkage.CVRods{i}(j+1).xi_star;
        gi      = Linkage.VLinks(Linkage.LinkIndex(i)).gi{j};

        Xs  = Linkage.CVRods{i}(j+1).Xs;
        nip = Linkage.CVRods{i}(j+1).nip;

        if Linkage.Z_order == 4
            Phi_Z1 = Linkage.CVRods{i}(j+1).Phi_Z1;
            Phi_Z2 = Linkage.CVRods{i}(j+1).Phi_Z2;
        else
            Phi_Z = Linkage.CVRods{i}(j+1).Phi_Z;
        end

        % enter element
        g_here = g_here * gi;

        if  s <= s_acc
            g = g_here;
            return
        end

        % -------------------------
        % Gauss/Zanna integration
        % -------------------------
        for ii = 2:nip

            H_full = (Xs(ii) - Xs(ii-1)) * Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};

            s_next = s_acc + H_full;

            % truncate step if s lies inside this interval
            if s <= s_next
                H = s - s_acc;
            else
                H = H_full;
            end   
            
            %smooth transition!
            % x = s - s_acc;
            % eps_s = 1e-6; % tune this
            % 
            % d = x - H_full;
            % H = 0.5*(x + H_full - sqrt(d*d + eps_s^2));


            % ---- compute strain ----
            if Linkage.Z_order == 4

                xi_Z1 = xi_star(6*(ii-2)+1:6*(ii-1),2);
                xi_Z2 = xi_star(6*(ii-2)+1:6*(ii-1),3);

                if dof_here > 0
                    Phi_Z1here = Phi_Z1(6*(ii-2)+1:6*(ii-1),:);
                    Phi_Z2here = Phi_Z2(6*(ii-2)+1:6*(ii-1),:);

                    xi_Z1 = Phi_Z1here*q_here + xi_Z1;
                    xi_Z2 = Phi_Z2here*q_here + xi_Z2;
                end

                ad_xi = dinamico_adj(xi_Z1);
                Omega = (H/2)*(xi_Z1 + xi_Z2) + ...
                        ((sqrt(3)*H^2)/12)*ad_xi*xi_Z2;

            else
                xi_Z = xi_star(6*(ii-2)+1:6*(ii-1),4);

                if dof_here > 0
                    Phi_Zhere = Phi_Z(6*(ii-2)+1:6*(ii-1),:);
                    xi_Z = Phi_Zhere*q_here + xi_Z;
                end

                Omega = H * xi_Z;
            end

            % exponential map
            gh = variable_expmap_g(Omega);
            g_here = g_here * gh;

            % -------------------------
            % STOP CONDITION
            % -------------------------
            if  s <= s_next
                g = g_here;
                return
            end

            s_acc = s_next;

        end

        gf = Linkage.VLinks(Linkage.LinkIndex(i)).gf{j};
        g_here = g_here * gf;

        if s <= s_acc
            g = g_here;
            return
        end

        dof_start = dof_start + dof_here;

    end

    % update tip
    g_tip((i-1)*4+1:i*4,:) = g_here;

end

% full result
g = g_here;

end