function Jd = JacobianDotAtS(Linkage, q, qd, s)

if isrow(q);  q  = q';  end
if isrow(qd); qd = qd'; end

N       = Linkage.N;
ndof    = Linkage.ndof;

g_ini   = Linkage.g_ini;
iLpre   = Linkage.iLpre;

g_tip   = repmat(eye(4),N,1);
Jd_tip  = zeros(6*N,ndof);
eta_tip = zeros(6*N,1);

Jd_here  = zeros(6,ndof);
eta_here = zeros(6,1);

dof_start = 1;

for i = 1:N

    % --- base transform ---
    if iLpre(i)>0
        g_here       = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:) * g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        Jd_here      = Ad_g_ini_inv * Jd_tip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        eta_here     = Ad_g_ini_inv * eta_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else
        g_here  = g_ini((i-1)*4+1:i*4,:);
        Jd_here = zeros(6,ndof);
        eta_here = zeros(6,1);
    end

    % =========================
    % JOINT (same as original)
    % =========================
    dof_here = Linkage.CVRods{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = qd(dof_start:dof_start+dof_here-1);

    Phi_here = Linkage.CVRods{i}(1).Phi;
    xi_star  = Linkage.CVRods{i}(1).xi_star;

    if dof_here > 0
        xi  = Phi_here*q_here + xi_star;
        xid = Phi_here*qd_here;

        [g_joint,Tg,Tgd] = variable_expmap_gTgTgd_mex(xi,xid);

        S_here = zeros(6,ndof);
        Sd_here = zeros(6,ndof);

        S_here(:,dof_start:dof_start+dof_here-1)  = Tg*Phi_here;
        Sd_here(:,dof_start:dof_start+dof_here-1) = ...
            dinamico_adj(eta_here)*Tg*Phi_here + Tgd*Phi_here;
    else
        g_joint = eye(4);
        S_here  = zeros(6,ndof);
        Sd_here = zeros(6,ndof);
    end

    g_here = g_here * g_joint;

    Ad = dinamico_Adjoint(ginv(g_joint));
    Jd_here = Ad*(Jd_here + Sd_here);
    eta_here = Ad*(eta_here + S_here(:,dof_start:dof_start+dof_here-1)*qd_here);

    dof_start = dof_start + dof_here;

    % =========================
    % SOFT LINK
    % =========================
    for j = 1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1

        ld  = Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};
        Xs  = Linkage.CVRods{i}(j+1).Xs;
        nip = Linkage.CVRods{i}(j+1).nip;

        % find segment where s lies
        ii_s = find(Xs >= s,1);
        if isempty(ii_s); ii_s = nip; end

        % --- propagate fully up to segment before s ---
        for ii = 2:ii_s-1
            H = (Xs(ii)-Xs(ii-1))*ld;

            [Jd_here, eta_here, g_here] = ...
                stepJacobianDot(Linkage,i,j,ii,H,q,qd,dof_start,...
                                Jd_here,eta_here,g_here);
        end

        % =========================
        % FINAL PARTIAL STEP
        % =========================
        if ii_s >= 2
            H = (s - Xs(ii_s-1)) * ld;

            [Jd_here, eta_here, g_here] = ...
                stepJacobianDot(Linkage,i,j,ii_s,H,q,qd,dof_start,...
                                Jd_here,eta_here,g_here);
        end

        Jd = Jd_here;
        return
    end

    g_tip((i-1)*4+1:i*4,:)   = g_here;
    Jd_tip((i-1)*6+1:i*6,:)  = Jd_here;
    eta_tip((i-1)*6+1:i*6,:) = eta_here;
end

end

function [Jd_here, eta_here, g_here] = stepJacobianDot(...
    Linkage,i,j,ii,H,q,qd,dof_start,Jd_here,eta_here,g_here)

ndof = Linkage.ndof;

dof_here = Linkage.CVRods{i}(j+1).dof;
q_here   = q(dof_start:dof_start+dof_here-1);
qd_here  = qd(dof_start:dof_start+dof_here-1);

xi_star = Linkage.CVRods{i}(j+1).xi_star;

Phi_Z = Linkage.CVRods{i}(j+1).Phi_Z(6*(ii-2)+1:6*(ii-1),:);

xi = xi_star(6*(ii-2)+1:6*(ii-1),4);

if dof_here > 0
    xi = Phi_Z*q_here + xi;
    Z  = H * Phi_Z;
    Omegad = Z * qd_here;
else
    Z = H * Phi_Z;
    Omegad = zeros(6,1);
end

Omega = H * xi;

[gh,Tg,Tgd] = variable_expmap_gTgTgd_mex(Omega,Omegad);

S  = zeros(6,ndof);
Sd = zeros(6,ndof);

S(:,dof_start:dof_start+dof_here-1)  = Tg*Z;
Sd(:,dof_start:dof_start+dof_here-1) = ...
    dinamico_adj(eta_here)*S(:,dof_start:dof_start+dof_here-1) + Tgd*Z;

g_here = g_here * gh;

Ad = dinamico_Adjoint(ginv(gh));

Jd_here  = Ad * (Jd_here + Sd);
eta_here = Ad * (eta_here + S(:,dof_start:dof_start+dof_here-1)*qd_here);

end