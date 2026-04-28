function J = JacobianAtS(Linkage,q,s)
% Jacobian evaluated up to arc-length position s ∈ [0,1]

if isrow(q)
    q = q';
end

ndof = Linkage.ndof;
N = Linkage.N;

g_ini = Linkage.g_ini;
g_tip = repmat(eye(4),N,1);
J_tip = zeros(6*N,ndof);
iLpre = Linkage.iLpre;

g_here = eye(4);
J_here = zeros(6,ndof);

dof_start = 1;
s_acc = 0;

J = zeros(6,ndof);
done = false;

for i = 1:N

    % -------------------------
    % base transform
    % -------------------------
    if iLpre(i)>0
        g_here = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:) * ...
                 g_ini((i-1)*4+1:i*4,:);
        Ad_prev = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        J_here = Ad_prev * J_tip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
    else
        g_here = g_ini((i-1)*4+1:i*4,:);
        J_here = zeros(6,ndof);
    end

    % -------------------------
    % joint
    % -------------------------
    dof_here = Linkage.CVRods{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    Phi_here = Linkage.CVRods{i}(1).Phi;
    xi_star  = Linkage.CVRods{i}(1).xi_star;

    if dof_here == 0
        g_joint = eye(4);
        S_here = zeros(6,ndof);
    else
        xi = Phi_here*q_here + xi_star;
        [g_joint,Tg] = variable_expmap_gTg(xi);

        S_here = zeros(6,ndof);
        S_here(:,dof_start:dof_start+dof_here-1) = Tg*Phi_here;
    end

    g_here = g_here * g_joint;
    J_here = dinamico_Adjoint(ginv(g_joint))*(J_here + S_here);

    % STOP check
    if s_acc >= s
        J = J_here;
        return
    end

    % rigid link
    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype == 'r'

        gi = Linkage.VLinks(Linkage.LinkIndex(i)).gi;
        gf = Linkage.VLinks(Linkage.LinkIndex(i)).gf;

        g_here = g_here * gi;
        J_here = dinamico_Adjoint(ginv(gi)) * J_here;

        H_rigid = Linkage.VLinks(Linkage.LinkIndex(i)).ld{1}; % if defined, else 1
        s_next = s_acc + H_rigid;

        if s <= s_next
            J = J_here;
            return
        end

        s_acc = s_next;

        g_here = g_here * gf;
        J_here = dinamico_Adjoint(ginv(gf)) * J_here;
    end

    dof_start = dof_start + dof_here;

    % -------------------------
    % flexible pieces
    % -------------------------
    for j = 1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1

        dof_here = Linkage.CVRods{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);

        xi_star = Linkage.CVRods{i}(j+1).xi_star;
        gi      = Linkage.VLinks(Linkage.LinkIndex(i)).gi{j};

        Xs  = Linkage.CVRods{i}(j+1).Xs;
        nip = Linkage.CVRods{i}(j+1).nip;
        ld  = Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};

        if Linkage.Z_order==4
            Phi_Z1 = Linkage.CVRods{i}(j+1).Phi_Z1;
            Phi_Z2 = Linkage.CVRods{i}(j+1).Phi_Z2;
        else
            Phi_Z = Linkage.CVRods{i}(j+1).Phi_Z;
        end

        g_here = g_here * gi;
        J_here = dinamico_Adjoint(ginv(gi)) * J_here;

        if s_acc >= s
            J = J_here;
            return
        end

        % -------------------------
        % Zanna integration
        % -------------------------
        for ii = 2:nip

            H_full = (Xs(ii)-Xs(ii-1)) * ld;
            s_next = s_acc + H_full;

            % truncate last step if needed
            if s <= s_next
                H = s - s_acc;
                done = true;
            else
                H = H_full;
            end

            % -------------------------
            % strain + sensitivity
            % -------------------------
            if Linkage.Z_order==4

                xi_Z1 = xi_star(6*(ii-2)+1:6*(ii-1),2);
                xi_Z2 = xi_star(6*(ii-2)+1:6*(ii-1),3);

                Phi_Z1here = Phi_Z1(6*(ii-2)+1:6*(ii-1),:);
                Phi_Z2here = Phi_Z2(6*(ii-2)+1:6*(ii-1),:);

                if dof_here>0
                    xi_Z1 = Phi_Z1here*q_here + xi_Z1;
                    xi_Z2 = Phi_Z2here*q_here + xi_Z2;
                end

                ad_xi = dinamico_adj(xi_Z1);

                Z_here = (H/2)*(Phi_Z1here+Phi_Z2here) + ...
                         ((sqrt(3)*H^2)/12)* ...
                         (ad_xi*Phi_Z2here - dinamico_adj(xi_Z2)*Phi_Z1here);

                Omega = (H/2)*(xi_Z1 + xi_Z2) + ...
                        ((sqrt(3)*H^2)/12)*ad_xi*xi_Z2;

            else

                xi_Z = xi_star(6*(ii-2)+1:6*(ii-1),4);
                Phi_Zhere = Phi_Z(6*(ii-2)+1:6*(ii-1),:);

                if dof_here>0
                    xi_Z = Phi_Zhere*q_here + xi_Z;
                end

                Z_here = H * Phi_Zhere;
                Omega  = H * xi_Z;
            end

            % exponential map
            [gh,Tg] = variable_expmap_gTg(Omega);

            S_here = zeros(6,ndof);
            S_here(:,dof_start:dof_start+dof_here-1) = Tg * Z_here;

            g_here = g_here * gh;
            J_here = dinamico_Adjoint(ginv(gh))*(J_here + S_here);

            % STOP exactly at s
            if done
                J = J_here;
                return
            end

            s_acc = s_next;

        end

        gf = Linkage.VLinks(Linkage.LinkIndex(i)).gf{j};
        g_here = g_here * gf;
        J_here = dinamico_Adjoint(ginv(gf)) * J_here;

        dof_start = dof_start + dof_here;

        if s_acc >= s
            J = J_here;
            return
        end

    end

    g_tip((i-1)*4+1:i*4,:) = g_here;
    J_tip((i-1)*6+1:i*6,:) = J_here;

end

J = J_here;

end