function g = InterpolateSE3(gk, gkp1, sk, skp1, s)

alpha = (s - sk)/(skp1 - sk);
alpha = max(0,min(1,alpha));

g_rel = gk \ gkp1;

xi_rel = variable_logmap_g(g_rel);

g = gk * variable_expmap_g(alpha * xi_rel);

end