    function err = FKErrorAt0(Xpos, T_0_fixed)
        %Function that describes error between floating frame coordinates
        %and fixed frame
        T_0 =variable_expmap_g(Xpos(1:6));

         dT = ginv(T_0) * T_0_fixed ;
         err = piecewise_logmap(dT);
    end