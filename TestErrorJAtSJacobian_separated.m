
%generate constraint satisfying carriage and rail configs
s_range = 0.1:1:5.3;
q_def_z = zeros(length(s_range),1);
q_b0 = 0.1*ones(RailLinkage.ndof,1);
for i = 1:length(s_range)
    g_s = FwdKinematicsAtS(RailLinkage,q_b0,s_range(i)); %SE(3)
    q_m0 = piecewise_logmap(g_s); % carriage coordinates (in tangent space of origin)

    %define the initial system state (position and vel)
    q0 = [q_b0;q_m0];

    TestErrorJAtSJacobian1(RailLinkage,RailCarriage,T_0_fixed,T_L_fixed,s_range(i),q0,zeros(ndof,1))
end
T_full_init = FwdKinematics(RailLinkage, x0(1:RailLinkage.ndof));
T_L_fixed = T_full_init(end-3:end, :);
T_0_fixed = T_full_init(1:4,:);

function TestErrorJAtSJacobian1(Linkage,Carriage,T_0_fixed,T_L_fixed,s,q,qd)
%disp(q)
    %disp(s)
[e,J] = ErrorJAtS_both_ends(Linkage,Carriage,s,q,qd);
   
% [err0,J0,Jd0,e,J,JdL] = ErrorDynamicsAt0andL( ...
%     Linkage, Carriage, q, qd, T_0_fixed, T_L_fixed);

nq = length(q);

eps_fd = 1e-7;

fprintf('\nTesting Jacobian:\n');

for i = 1:nq

    dq = zeros(nq,1);
    dq(i) = eps_fd;

    ep = ErrorJAtS_both_ends(Linkage,Carriage,s,q+dq,qd);

    % ep = ErrorDynamicsAt0andL( ...
    %     Linkage, Carriage, q+dq, qd, T_0_fixed, T_L_fixed);

    em = ErrorJAtS_both_ends(Linkage,Carriage,s,q-dq,qd);

    % em = ErrorDynamicsAt0andL( ...
    %     Linkage, Carriage, q-dq, qd, T_0_fixed, T_L_fixed);

    J_fd_col = (ep-em)/(2*eps_fd);

    err = norm(J_fd_col - J(:,i));

    fprintf('%3d  abs err = %.3e\n',i,err);

end

end


L = RailCarriage.VLinks(1).L;

T_1_offset = [eye(3) [-L/2;0;0];0 0 0 1];
T_2_offset = [eye(3) [L/2;0;0];0 0 0 1];

q_test = ones(6,1);
[T,Jm] = variable_expmap_gTgTgd(q_test,zeros(6,1));

for i = 1:6

    dq = zeros(6,1);
    dq(i) = 1e-7;

    Tp = variable_expmap_g(q_test+dq)*T_1_offset;
    Tm = variable_expmap_g(q_test-dq)*T_1_offset;

    delta = piecewise_logmap(ginv(Tm)*Tp);

    Jfd(:,i) = delta/(2e-7);

end
Jfd