<<<<<<< HEAD
function u = PD_control(Linkage, s, q, qd, controller)

qd_mass = qd(Linkage.ndof+1:end);

u = controller.kp*(controller.s_cmd - s) +  controller.kd*(controller.sd_cmd - qd_mass(4));
=======
function u = PD_carriage_control(Linkage, s, q, qd, PD_controller)

qd_mass = qd(Linkage.ndof+1:end);

u = PD_controller.kp*(PD_controller.s_cmd - s) + PD_controller.kd*(PD_controller.sd_cmd - qd_mass(4));
>>>>>>> cc57c9f2fb11373d2b25f7d18c94213a2451e0e9

end