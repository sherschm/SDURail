function u = PD_control(Linkage, s, q, qd, controller)

qd_mass = qd(Linkage.ndof+1:end);

u = controller.kp*(controller.s_cmd - s) +  controller.kd*(controller.sd_cmd - qd_mass(4));

end