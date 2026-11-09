
q_b = x_out(end,1:RailLinkage.ndof)';
qd_b = ones(RailLinkage.ndof,1);

[err0,J0] = ErrorDynamicsAt0andL( ...
    RailLinkage,RailCarriage,q_b,qd_b,T_0_fixed,T_L_fixed);

eps_fd = 1e-7;

for i = 1:length(q_b)

    dq = zeros(length(q_b),1);
    dq(i) = eps_fd;

    [err0p,~,~,~,~,~] = ErrorDynamicsAt0andL( ...
        RailLinkage,RailCarriage,...
        q_b+dq,qd_b,...
        T_0_fixed,T_L_fixed);

    [err0m,~,~,~,~,~] = ErrorDynamicsAt0andL( ...
        RailLinkage,RailCarriage,...
        q_b-dq,qd_b,...
        T_0_fixed,T_L_fixed);

    Jfd_col = (err0p - err0m)/(2*eps_fd);

    err = norm(Jfd_col - J0(:,i));

    fprintf('q(%d): %.3e\n',i,err);

end

dq = randn(length(q_b),1);
dq = dq/norm(dq);

h = 1e-7;

[err0p,~,~,~,~,~] = ErrorDynamicsAt0andL( ...
    RailLinkage,RailCarriage,...
    q_b+h*dq,...
    qd_b,...
    T_0_fixed,T_L_fixed);

[err0m,~,~,~,~,~] = ErrorDynamicsAt0andL( ...
    RailLinkage,RailCarriage,...
    q_b-h*dq,...
    qd_b,...
    T_0_fixed,T_L_fixed);

de_fd = (err0p-err0m)/(2*h);

de_J = J0(:,1:length(q_b))*dq;

fprintf('Error = %.3e\n', ...
    norm(de_fd-de_J));


[err0,J0,~,errL,JL,~] = ErrorDynamicsAt0andL( ...
    RailLinkage,RailCarriage,...
    q_b,qd_b,...
    T_0_fixed,T_L_fixed);

Janalytic = [J0(:,1:length(q_b));
    JL(:,1:length(q_b))];

Jfd = NumericalJacobian_ErrorDynamicsAt0andL( ...
    RailLinkage,RailCarriage,...
    q_b,qd_b,...
    T_0_fixed,T_L_fixed);

fprintf('||Jfd-J|| = %.3e\n', norm(Jfd-Janalytic));
fprintf('relative = %.3e\n', ...
    norm(Jfd-Janalytic)/max(norm(Jfd),1e-12));

function Jfd = NumericalJacobian_ErrorDynamicsAt0andL( ...
    Linkage,Carriage,q_b,qd_b,T_0_fixed,T_L_fixed)

eps_fd = 1e-4;

[err0,~,~,errL,~,~] = ErrorDynamicsAt0andL( ...
    Linkage,Carriage,...
    q_b,qd_b,...
    T_0_fixed,T_L_fixed);

e0 = [err0;
    errL];

m = length(e0);
n = length(q_b);

Jfd = zeros(m,n);

for i = 1:n

    dq = zeros(n,1);
    dq(i) = eps_fd;

    [err0p,~,~,errLp,~,~] = ErrorDynamicsAt0andL( ...
        Linkage,Carriage,...
        q_b+dq,qd_b,...
        T_0_fixed,T_L_fixed);

    [err0m,~,~,errLm,~,~] = ErrorDynamicsAt0andL( ...
        Linkage,Carriage,...
        q_b-dq,qd_b,...
        T_0_fixed,T_L_fixed);

    ep = [err0p;
        errLp];

    em = [err0m;
        errLm];

    Jfd(:,i) = (ep-em)/(2*eps_fd);

end
end