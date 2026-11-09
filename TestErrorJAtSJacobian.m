
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
[e,J] = ErrorJAtS(Linkage,Carriage,s,q,qd);
   
% [err0,J0,Jd0,e,J,JdL] = ErrorDynamicsAt0andL( ...
%     Linkage, Carriage, q, qd, T_0_fixed, T_L_fixed);

nq = length(q);

eps_fd = 1e-7;

fprintf('\nTesting Jacobian:\n');

for i = 1:nq

    dq = zeros(nq,1);
    dq(i) = eps_fd;

    ep = ErrorJAtS(Linkage,Carriage,s,q+dq,qd);

    % ep = ErrorDynamicsAt0andL( ...
    %     Linkage, Carriage, q+dq, qd, T_0_fixed, T_L_fixed);

    em = ErrorJAtS(Linkage,Carriage,s,q-dq,qd);

    % em = ErrorDynamicsAt0andL( ...
    %     Linkage, Carriage, q-dq, qd, T_0_fixed, T_L_fixed);

    J_fd_col = (ep-em)/(2*eps_fd);

    err = norm(J_fd_col - J(:,i));

    fprintf('%3d  abs err = %.3e\n',i,err);

end

end

q_test = ones(6,1);
[T,Jm] = variable_expmap_gTgTgd(q_test,zeros(6,1));

for i = 1:6

    dq = zeros(6,1);
    dq(i) = 1e-7;

    Tp = variable_expmap_g(q_test+dq);
    Tm = variable_expmap_g(q_test-dq);

    delta = piecewise_logmap(ginv(Tm)*Tp);

    Jfd(:,i) = delta/(2e-7);

end
Jfd


qd_test = ones(6,1);

[T,Jm,Jmd_wrong] = variable_expmap_gTgTgd(q_test,qd_test);

h = 1e-7;

qp = q_test + h*qd_test;
qm = q_test - h*qd_test;

[~,Jp_wrong] = variable_expmap_gTgTgd(qp,qd_test);
[~,Jm_wrong] = variable_expmap_gTgTgd(qm,qd_test);

% Apply your correction
Jp = zeros(6,6);
Jp(1:3,1:3) = Jp_wrong(1:3,1:3)';
Jp(4:6,4:6) = Jp_wrong(4:6,4:6)';
Jp(4:6,1:3) = Jp_wrong(4:6,1:3)';

Jm_fd = zeros(6,6);
Jm_fd(1:3,1:3) = Jm_wrong(1:3,1:3)';
Jm_fd(4:6,4:6) = Jm_wrong(4:6,4:6)';
Jm_fd(4:6,1:3) = Jm_wrong(4:6,1:3)';

Jmd_fd = (Jp - Jm_fd)/(2*h);

Jmd = zeros(6,6);
Jmd(1:3,1:3) = Jmd_wrong(1:3,1:3)';
Jmd(4:6,4:6) = Jmd_wrong(4:6,4:6)';
Jmd(4:6,1:3) = Jmd_wrong(4:6,1:3)';

disp(norm(Jmd_fd - Jmd))