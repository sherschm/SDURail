function [ds_dqb, ds_dqm] = ProjectSJacobianFD(Linkage,q_b,q_mass,s)

eps_fd = 1e-4;

nqb = length(q_b);
nqm = length(q_mass);

ds_dqb = zeros(1,nqb);
ds_dqm = zeros(1,nqm);

% -------------------------
% ds/dq_b
% -------------------------
for i = 1:nqb

    dq = zeros(nqb,1);
    dq(i) = eps_fd;

    sp = ProjectS(Linkage,q_b+dq,q_mass,s);
    sm = ProjectS(Linkage,q_b-dq,q_mass,s);

    ds_dqb(i) = (sp-sm)/(2*eps_fd);

end

% -------------------------
% ds/dq_mass
% -------------------------
for i = 1:nqm

    dq = zeros(nqm,1);
    dq(i) = eps_fd;

    sp = ProjectS(Linkage,q_b,q_mass+dq,s);
    sm = ProjectS(Linkage,q_b,q_mass-dq,s);

    ds_dqm(i) = (sp-sm)/(2*eps_fd);

end

end