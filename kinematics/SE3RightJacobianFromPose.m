function Jr = SE3RightJacobianFromPose(xi)

[g_joint,Jr] = variable_expmap_gTg(xi);
  
%Jr = dinamico_Adjoint(ginv(g_joint))* Tg;

end