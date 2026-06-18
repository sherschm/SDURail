
function Jd = ErrorJdAtS_FD(Linkage, Carriage, s, q, qd) 
        
    eps_fd = 1e-6;

    q_plus  = q + eps_fd*qd;
    q_minus = q - eps_fd*qd;
    
    [~, J_plus]  = ErrorJAtS(Linkage,Carriage,s,q_plus);
    [~, J_minus] = ErrorJAtS(Linkage,Carriage,s,q_minus);
    
    Jd = (J_plus - J_minus)/(2*eps_fd);

end
