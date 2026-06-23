function u = PD_carriage_control(Linkage, s, q, qd, PD_controller)

qd_mass = qd(Linkage.ndof+1:end);

u = PD_controller.kp*(PD_controller.s_cmd - s) + PD_controller.kd*(PD_controller.sd_cmd - qd_mass(4));

end