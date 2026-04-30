function err = carriage_constraint_err(Linkage, s, q)
        % Split generalized coordinates
        n_beam = Linkage.ndof;
        n_mass = 6;
        
        q_b = q(1:n_beam);
        q_mass = q(n_beam+1:n_beam+n_mass);
        
        err_full = piecewise_logmap(FwdKinematicsAtS(Linkage,q_b,s))-q_mass;

        
        if Carriage.fixed == true
            err = err_full;
        else
            err = [err_full(1:3);err_full(5:6)];
        end
    end