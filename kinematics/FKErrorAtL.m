    function err = FKErrorAtL(Linkage, Xpos, T_L_fixed, n_beam)
        %Function that describes error between end pose coordinates
        %and fixed end position.

        q_b = Xpos(1:n_beam);
        
        T_L_full = FwdKinematics(Linkage, q_b);
        T_L = T_L_full(end-3:end,:);
        

         dT = ginv(T_L) * T_L_fixed ;
         err = piecewise_logmap(dT);

    end