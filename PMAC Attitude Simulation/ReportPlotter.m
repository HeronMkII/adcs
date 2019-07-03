
form = 'b';    % changes line format

%% omega and Euler angles
f1 = figure(1); set(f1, 'defaulttextinterpreter','latex')

subplot(3,2,1); hold on;
plot(t_axis,X_qua(1,:),form);
title('Euler Param 1 vs. time')
xlabel('Orbit'); ylabel('$\epsilon1$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,3); hold on;
plot(t_axis,X_qua(2,:),form);
title('Euler Param 2 vs. time')
xlabel('Orbit'); ylabel('$\epsilon2$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,5); hold on;
plot(t_axis,X_qua(3,:),form);
title('Euler Param 3 vs. time')
xlabel('Orbit'); ylabel('$\epsilon3$');
grid on; grid minor; xlim([0 3]);

subplot(3,2,2) ; hold on;
plot(t_axis,X_ang(1,:),form);
title('Angular Velocity 1 vs. time')
xlabel('Orbit'); ylabel('$\omega_1$ (rad/s)');
grid on; grid minor; xlim([0 3]);
subplot(3,2,4) ; hold on;
plot(t_axis,X_ang(2,:),form);
title('Angular Velocity 2 vs. time')
xlabel('Orbit'); ylabel('$\omega_2$ (rad/s)');
grid on; grid minor; xlim([0 3]);
subplot(3,2,6) ; hold on;
plot(t_axis,X_ang(3,:),form);
title('Angular Velocity 3 vs. time')
xlabel('Orbit'); ylabel('$\omega_3$ (rad/s)');
grid on; grid minor; xlim([0 3]);


%% tau_mag and tau_imp 
f2 = figure(2); set(f2, 'DefaultTextInterpreter','latex');

subplot(3,2,1); hold on;
plot(t_axis,u_ct_ary(1,:),form);
title('Cont input (dipole) $u_{ct,1}$ 1 vs. time')
xlabel('Orbit'); ylabel('$m_1$ $(Am^2)$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,3); hold on;
plot(t_axis,u_ct_ary(2,:),form);
title('Cont input (dipole) $u_{ct,2}$ vs. time')
xlabel('Orbit'); ylabel('$m_2$ $(Am^2)$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,5); hold on;
plot(t_axis,u_ct_ary(3,:),form);
title('Cont input (dipole) $u_{ct,3}$ vs. time')
xlabel('Orbit'); ylabel('$m_3$ $(Am^2)$');
grid on; grid minor; xlim([0 3]);

subplot(3,2,2); hold on;
plot(t_axis,u_ds_ary(1,:),form);
title('Discrete input $u_{ds,1}$ vs. time')
xlabel('Orbit'); ylabel('$\tau_{imp,1}$ $(Nm)$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,4); hold on;
plot(t_axis,u_ds_ary(2,:),form);
title('Discrete input $u_{ds,2}$ vs. time')
xlabel('Orbit'); ylabel('$\tau_{imp,2}$ $(Nm)$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,6); hold on;
plot(t_axis,u_ds_ary(3,:),form);
title('Discrete input $u_{ds,3}$ vs. time')
xlabel('Orbit'); ylabel('$\tau_{imp,3}$ $(Nm)$');
grid on; grid minor; xlim([0 3]);

%% Energy and wheel speed
f3 = figure(3); set(f3, 'DefaultTextInterpreter','latex');
subplot(3,2,1); hold on;
plot(t_axis,MagTorq.cur(1,:),form);
title('MagTorquer Current 1 vs. time')
xlabel('Orbit'); ylabel('$I_1$ $(A)$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,3); hold on; 
plot(t_axis,MagTorq.cur(2,:),form);
title('MagTorquer Current 2 vs. time')
xlabel('Orbit'); ylabel('$I_2$ $(A)$');
grid on; grid minor; xlim([0 3]);
subplot(3,2,5); hold on;
plot(t_axis,MagTorq.cur(3,:),form);
title('MagTorquer Current 3 vs. time')
xlabel('Orbit'); ylabel('$I_3$ $(A)$');
grid on; grid minor; xlim([0 3]);

subplot(3,2,2); plot(t_axis,Wheel.AngSpd(1,:)*60/2/pi,form); hold on;
title('Reaction Wheel 1 vs. time')
xlabel('Orbit'); ylabel('RPM');
grid on; grid minor; xlim([0 3]);
subplot(3,2,4); plot(t_axis,Wheel.AngSpd(2,:)*60/2/pi,form); hold on;
title('Reaction Wheel 2 vs. time')
xlabel('Orbit'); ylabel('RPM');
grid on; grid minor; xlim([0 3]);
subplot(3,2,6); plot(t_axis,Wheel.AngSpd(3,:)*60/2/pi,form); hold on;
title('Reaction Wheel 3 vs. time')
xlabel('Orbit'); ylabel('RPM');
grid on; grid minor; xlim([0 3]);


%% u_mag and u_ct overlapped (individual)
f4 = figure; set(f4,'DefaultTextInterpreter','latex');

subplot(3,1,1)
plot(t_axis,torq.mag(1,:)); ylabel('$\tau_{mag,1}$ $(N\cdot m)$'); ylim(4e-6*[-1 1]);
yyaxis right; 
plot(t_axis,torq.imp(1,:)); ylabel('$\tau_{imp,1}$ $(N\cdot m)$'); ylim(5e-5*[-1 1]);
xlabel('Orbit'); xlim([0 3]);
title('Controller efforts around axis 1'); grid on; grid minor;
subplot(3,1,2)
plot(t_axis,torq.mag(2,:)); ylabel('$\tau_{mag,2}$ $(N\cdot m)$'); ylim(4e-6*[-1 1]);
yyaxis right; 
plot(t_axis,torq.imp(2,:)); ylabel('$\tau_{imp,2}$ $(N\cdot m)$'); ylim(5e-5*[-1 1]);
xlabel('Orbit'); xlim([0 3]);
title('Controller efforts around axis 2'); grid on; grid minor;
subplot(3,1,3)
plot(t_axis,torq.mag(3,:)); ylabel('$\tau_{mag,3}$ $(N\cdot m)$'); ylim(4e-6*[-1 1]);
yyaxis right; 
plot(t_axis,torq.imp(3,:)); ylabel('$\tau_{imp,3}$ $(N\cdot m)$'); ylim(5e-5*[-1 1]);
xlabel('Orbit'); xlim([0 3]);
title('Controller efforts around axis 3'); grid on; grid minor;
 