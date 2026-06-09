%Function that calculates the jacobian at every significant points
% in SI units (Last modified 18.04.2023, Anup)
function J = JacobianAtPoint(Linkage,q,i_here,j_here)%i_here is link index, j_here is division (0 for joint)

if isrow(q)
    q=q';
end

full = false;
if nargin==4
    if j_here==0
        nsig = 1;
    else
        nsig = Linkage.CVRods{i_here}(j_here+1).nip; %j_here>1 is allowed only for soft links
    end
elseif nargin==3
    if Linkage.VLinks(Linkage.LinkIndex(i_here)).linktype=='s'
        nsig = 1;
        for j=1:Linkage.VLinks(Linkage.LinkIndex(i_here)).npie-1
            nsig = nsig+Linkage.CVRods{i_here}(j+1).nip;
        end
    else
        nsig = 2; %joint and CM
    end
else
    full  = true;
    nsig  = Linkage.nsig;
end
ndof  = Linkage.ndof;
J     = zeros(6*nsig,ndof);
i_sig = 1;

N         = Linkage.N;

g_ini = Linkage.g_ini;
g_tip = repmat(eye(4),N,1);
J_tip = zeros(6*N,ndof);
iLpre = Linkage.iLpre;

dof_start = 1; %starting dof of current piece

for i = 1:N
    
    if iLpre(i)>0
        g_here       = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        J_here       = Ad_g_ini_inv*J_tip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
    else
        g_here = g_ini((i-1)*4+1:i*4,:);
        J_here = zeros(6,ndof);
    end
    
    %Joint
    dof_here = Linkage.CVRods{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    Phi_here = Linkage.CVRods{i}(1).Phi;
    xi_star  = Linkage.CVRods{i}(1).xi_star;
    
    if dof_here == 0 %fixed joint (N)
        g_joint  = eye(4);
        S_here = zeros(6,ndof);
    else
        xi           = Phi_here*q_here+xi_star;
        [g_joint,Tg] = variable_expmap_gTg(xi);

        S_here                                   = zeros(6,ndof);
        S_here(:,dof_start:dof_start+dof_here-1) = Tg*Phi_here;
    end
    
    %updating g, Jacobian, Jacobian_dot and eta
    g_here = g_here*g_joint;
    J_here = dinamico_Adjoint(ginv(g_joint))*(J_here+S_here);
    
    if full||(i==i_here&&nargin==3)||(i==i_here&&j_here==0)
        J((i_sig-1)*6+1:i_sig*6,:) = J_here;
        i_sig                      = i_sig+1;
    end
    
    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype == 'r'
        
        gi     = Linkage.VLinks(Linkage.LinkIndex(i)).gi;
        g_here = g_here*gi;
        J_here = dinamico_Adjoint(ginv(gi))*J_here;
        
        if full||(i==i_here&&nargin==3)
            J((i_sig-1)*6+1:i_sig*6,:) = J_here;
            i_sig                      = i_sig+1;
        end
        
        % bringing all quantities to the end of rigid link
        gf     = Linkage.VLinks(Linkage.LinkIndex(i)).gf;
        g_here = g_here*gf;
        J_here = dinamico_Adjoint(ginv(gf))*J_here;
    end
    J = J_here(1:6,1:6);
end
end

