function Tg_dot = variable_expmap_gTg_dot(Omega, Omega_dot)

eps = 1e-8;

[~,Tg0] = variable_expmap_gTg(Omega);
[~,Tg1] = variable_expmap_gTg(Omega + eps*Omega_dot);

Tg_dot = (Tg1 - Tg0)/eps;

end