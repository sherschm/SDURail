function u = PD_control(Linkage, s, q, qd, s_cmd, sd_cmd, k_p, k_d)

qd_mass = qd(Linkage.ndof+1:end);

u = k_p*(s_cmd - s) + k_d*(sd_cmd - qd_mass(4));

end