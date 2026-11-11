function [e,J,Jd] = ErrorJAtS_both_ends(Linkage, Carriage, s, q, qd, opts) 

arguments
    Linkage
    Carriage
    s
    q
    qd
    opts.fixed logical = false
end

if opts.fixed
    %S = eye(6);
    S =[1 0 0 0 0 0; ...
        0 0 0 1 0 0; ...
        0 0 0 0 1 0; ...
        0 0 0 0 0 1];
else
    % S = [1 0 0 0 0 0; ...
    %     0 1 0 0 0 0; ...
    %     0 0 1 0 0 0; ...
    %     0 0 0 0 1 0; ...
    %     0 0 0 0 0 1];

    S = [1 0 0 0 0 0; ...
        0 0 0 0 1 0; ...
        0 0 0 0 0 1];

end

%S= eye(6);

L = Carriage.VLinks(1).L;

% T_1_offset = [eye(3) [-L/2;0;0];0 0 0 1];
% T_2_offset = [eye(3) [L/2;0;0];0 0 0 1];
T_1_offset = eye(4);
%T_2_offset = ginv(ginv(Carriage.VLinks(1).gi)*Carriage.VLinks(1).gf);
T_2_offset = ginv(Carriage.VLinks(1).gf)*Carriage.VLinks(1).gi;

n_beam = Linkage.ndof;
n_mass = Carriage.ndof;
q_b = q(1:n_beam);
q_mass = q(n_beam+1:n_beam+n_mass);
qd_b = qd(1:n_beam);
qd_mass = qd(n_beam+1:n_beam+n_mass);

% --- Forward kinematics 
s1 = s-L/2;
s2 = s+L/2;

% s1 = s;
% s2 = s+L;

T_s1 = FwdKinematicsAtS(Linkage, q_b, s1);
T_s2 = FwdKinematicsAtS(Linkage, q_b, s2);

%Kinematics of the centre of the body:
[T_m,Jm_WRONG,Jmd_WRONG] = variable_expmap_gTgTgd(q_mass,qd_mass);

% %Kinematics of each ends
Tm_1 = T_m*T_1_offset;
Tm_2 = T_m*T_2_offset;

% Ad_offset_1 = dinamico_Adjoint(ginv(T_1_offset));
% 
% Jm_1_WRONG = Ad_offset_1*Jm_WRONG;
% Jmd_1_WRONG = Ad_offset_1*Jmd_WRONG;
% 
% Jm_1 = zeros(6,6);
% Jm_1(1:3,1:3) = Jm_1_WRONG(1:3,1:3)';
% Jm_1(4:6,4:6) = Jm_1_WRONG(4:6,4:6)';
% Jm_1(4:6,1:3) = Jm_1_WRONG(4:6,1:3)';
% Jmd_1 = zeros(6,6);
% Jmd_1(1:3,1:3) = Jmd_1_WRONG(1:3,1:3)';
% Jmd_1(4:6,4:6) = Jmd_1_WRONG(4:6,4:6)';
% Jmd_1(4:6,1:3) = Jmd_1_WRONG(4:6,1:3)';
% 

% Ad_offset_2 = dinamico_Adjoint(ginv(T_2_offset));
% Jm_2_WRONG = Ad_offset_2*Jm_WRONG;
% Jmd_2_WRONG = Ad_offset_2*Jmd_WRONG;
% 
% Jm_2 = zeros(6,6);
% Jm_2(1:3,1:3) = Jm_2_WRONG(1:3,1:3)';
% Jm_2(4:6,4:6) = Jm_2_WRONG(4:6,4:6)';
% Jm_2(4:6,1:3) = Jm_2_WRONG(4:6,1:3)';
% Jmd_2 = zeros(6,6);
% Jmd_2(1:3,1:3) = Jmd_2_WRONG(1:3,1:3)';
% Jmd_2(4:6,4:6) = Jmd_2_WRONG(4:6,4:6)';
% Jmd_2(4:6,1:3) = Jmd_2_WRONG(4:6,1:3)';

Ad_offset_1 = dinamico_Adjoint(ginv(T_1_offset));
Ad_offset_2 = dinamico_Adjoint(ginv(T_2_offset));
Jm = zeros(6,6);
Jm(1:3,1:3) = Jm_WRONG(1:3,1:3)';
Jm(4:6,4:6) = Jm_WRONG(4:6,4:6)';
Jm(4:6,1:3) = Jm_WRONG(4:6,1:3)';
Jmd = zeros(6,6);
Jmd(1:3,1:3) = Jmd_WRONG(1:3,1:3)';
Jmd(4:6,4:6) = Jmd_WRONG(4:6,4:6)';
Jmd(4:6,1:3) = Jmd_WRONG(4:6,1:3)';
Jm_1 = Ad_offset_1*Jm;
Jmd_1 = Ad_offset_1*Jmd;
Jm_2 = Ad_offset_2*Jm;
Jmd_2 = Ad_offset_2*Jmd;

g_sm_1 = ginv(T_s1) * Tm_1;
g_sm_2 = ginv(T_s2) * Tm_2;

% --- Jacobians in their natural frames
Js_1 = JacobianAtS(Linkage, q_b, s1);
Jsd_1 = JacobiandotAtS(Linkage, q_b,qd_b, s1);
Js_2 = JacobianAtS(Linkage, q_b, s2);
Jsd_2 = JacobiandotAtS(Linkage, q_b,qd_b, s2);

e  = [S*piecewise_logmap(g_sm_1);S*piecewise_logmap(g_sm_2)];
J =  [S*[-Js_1 Jm_1];S*[-Js_2 Jm_2]];
Jd =  [S*[-Jsd_1 Jmd_1];S*[-Jsd_2 Jmd_2]];
end
