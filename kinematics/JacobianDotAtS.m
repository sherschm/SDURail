function Jd_out = JacobianDotAtS(Linkage, q, qd, s)

if isrow(q);  q  = q';  end
if isrow(qd); qd = qd'; end

N     = Linkage.N;
ndof  = Linkage.ndof;
g_ini   = Linkage.g_ini;
iLpre   = Linkage.iLpre;
g_tip   = repmat(eye(4), N, 1);
Jd_tip  = zeros(6*N, ndof);
eta_tip = zeros(6*N, 1);

full  = true;
nsig  = Linkage.nsig;

Jd_out = zeros(6,ndof);

% Output (single point)
Jd = zeros(6, ndof);
i_sig     = 1;
dof_start = 1;
s_acc = 0;
done = false;


for i = 1:N
    
    if iLpre(i)>0
        g_here       = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        Jd_here      = Ad_g_ini_inv*Jd_tip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        eta_here     = Ad_g_ini_inv*eta_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        Jd_here  = zeros(6,ndof);
        eta_here   = zeros(6,1);
    end
    
    % ----------------------------
    % joint
    % ----------------------------
    dof_here = Linkage.CVRods{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = qd(dof_start:dof_start+dof_here-1);
    
    Phi_here = Linkage.CVRods{i}(1).Phi;
    xi_star  = Linkage.CVRods{i}(1).xi_star;
    
    if dof_here == 0
        g_joint = eye(4);
        S_here  = zeros(6,ndof);
        Sd_here = zeros(6,ndof);
    else
        xi               = Phi_here*q_here+xi_star;
        xid              = Phi_here*qd_here;
        [g_joint,Tg,Tgd] = variable_expmap_gTgTgd_mex(xi,xid);

        S_here                                    = zeros(6,ndof);
        S_here(:,dof_start:dof_start+dof_here-1)  = Tg*Phi_here;
        Sd_here                                   = zeros(6,ndof);
        Sd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*Tg*Phi_here+Tgd*Phi_here;
    end
    
    % propagate
    g_here = g_here * g_joint;
    Ad_inv = dinamico_Adjoint(ginv(g_joint));
    
    Jd_here  = Ad_inv * (Jd_here + Sd_here);
    eta_here = Ad_inv * (eta_here + ...
                S_here(:,dof_start:dof_start+dof_here-1)*qd_here);
    
    dof_start = dof_start + dof_here;
    
    if full||(i==i_here&&nargin==4)||(i==i_here&&j_here==0)
       % Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
        i_sig                       = i_sig+1;
    end
    % STOP check
    if s_acc >= s
        Jd = Jd_here;
        return
    end

    % ----------------------------
    % ONLY SOFT LINK HANDLED FOR s
    % ----------------------------
    for j = 1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
        
        dof_here = Linkage.CVRods{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        qd_here  = qd(dof_start:dof_start+dof_here-1);
        xi_star  = Linkage.CVRods{i}(j+1).xi_star;
        gi       = Linkage.VLinks(Linkage.LinkIndex(i)).gi{j};
        ld       = Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};
        Xs       = Linkage.CVRods{i}(j+1).Xs;
        nip      = Linkage.CVRods{i}(j+1).nip;
     
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here   = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        Jd_here   = Ad_gi_inv*Jd_here;
        eta_here  = Ad_gi_inv*eta_here;
        
        if full||(i==i_here&&nargin==4)||(i==i_here&&j==j_here)
            Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
            i_sig                       = i_sig+1;
        end 

        if s_acc >= s
            J = J_here;
            return
        end
        % ----------------------------
        % find interval containing s
        % ----------------------------
        for ii = 2:nip
            
            %H    = (Xs(ii)-Xs(ii-1))*ld;
            H_full = (Xs(ii)-Xs(ii-1)) * ld;
            s_next = s_acc + H_full;

            % truncate last step if needed
            if s <= s_next
                H = s - s_acc;
                done = true;
            else
                H = H_full;
            end

            
            if Linkage.Z_order==4
                
                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2); 
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                    
                Phi_Z1here  = Linkage.CVRods{i}(j+1).Phi_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                Phi_Z2here  = Linkage.CVRods{i}(j+1).Phi_Z2(6*(ii-2)+1:6*(ii-1),:);

                if dof_here>0
                    xi_Z1here = Phi_Z1here*q_here+xi_Z1here;
                    xi_Z2here = Phi_Z2here*q_here+xi_Z2here;

                    xid_Z1here  = Phi_Z1here*qd_here;
                    xid_Z2here  = Phi_Z2here*qd_here;

                    ad_xi_Z1here = dinamico_adj(xi_Z1here);

                    Z_here  = (H/2)*(Phi_Z1here+Phi_Z2here)+...
                                       ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*Phi_Z2here-dinamico_adj(xi_Z2here)*Phi_Z1here);

                    Zd_here = ((sqrt(3)*H^2)/12)*(dinamico_adj(xid_Z1here)*Phi_Z2here-dinamico_adj(xid_Z2here)*Phi_Z1here); 

                    Omegad_here   = Z_here*qd_here;
                else
                    ad_xi_Z1here = dinamico_adj(xi_Z1here);
                    Z_here  = (H/2)*(Phi_Z1here+Phi_Z2here)+...
                                       ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*Phi_Z2here-dinamico_adj(xi_Z2here)*Phi_Z1here);
                    Zd_here = zeros(6,dof_here); 
                    Omegad_here = zeros(6,1); 
                end

                Omega_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                             ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
                      
            else % order 2
                
                xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1),4);
                    
                Phi_Zhere  = Linkage.CVRods{i}(j+1).Phi_Z(6*(ii-2)+1:6*(ii-1),:);%note this step

                if dof_here>0
                    xi_Zhere = Phi_Zhere*q_here+xi_Zhere;

                    Z_here = H*Phi_Zhere;

                    Omegad_here     = Z_here*qd_here;
                else
                    Z_here = H*Phi_Zhere;
                    Omegad_here   = zeros(6,1); 
                end

                
                Omega_here  = H*xi_Zhere;

            end
            
            [gh,Tg,Tgd] = variable_expmap_gTgTgd_mex(Omega_here,Omegad_here); % mex code, C program
            
            S_here                                    = zeros(6,ndof);
            S_here(:,dof_start:dof_start+dof_here-1)  = Tg*Z_here;
            Sd_here                                   = zeros(6,ndof);
            Sd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*S_here(:,dof_start:dof_start+dof_here-1)+Tgd*Z_here;
            
            if Linkage.Z_order==4
                Sd_here(:,dof_start:dof_start+dof_here-1) = Sd_here(:,dof_start:dof_start+dof_here-1)+Tg*Zd_here;
            end
            
            %updating g, Jacobian, Jacobian_dot and eta
            g_here    = g_here*gh;
            Ad_gh_inv = dinamico_Adjoint(ginv(gh));
            Jd_here   = Ad_gh_inv*(Jd_here+Sd_here); %full
            eta_here  = Ad_gh_inv*(eta_here+S_here(:,dof_start:dof_start+dof_here-1)*qd_here);
            
            
            if full||(i==i_here&&nargin==4)||(i==i_here&&j==j_here)
                Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
                i_sig                       = i_sig+1;
            end
            

            % STOP exactly at s
            if done
                %J = J_here;
                Jd_out = Jd_here;
                return
            end
        
            % if s> Linkage.VLinks.L
            %     Jd_out = zeros(6,Linkage.ndof);%Jacobiandot()
            %     return
            % end
            s_acc = s_next;

        end
    
    g_tip((i-1)*4+1:i*4,:)   = g_here;
    Jd_tip((i-1)*6+1:i*6,:)  = Jd_here;
    eta_tip((i-1)*6+1:i*6,:) = eta_here;
end

end