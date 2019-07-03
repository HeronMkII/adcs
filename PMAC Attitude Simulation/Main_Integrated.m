clear;clc;close all;

%% ---Earths constants & Orbital Elements--- %%
% Earth Constants
Eth.G = 6.67e-11;                   % Gravitational constant
Eth.M = 5.972e24;                   % [kg] Mass of Earth
Eth.mu = Eth.G * Eth.M;             %[m3s-2] Gratitational constant
Eth.r_E = 6.3781e6;                 %[m] Earth's radius
    
% ISS orbit
orb.i = 51.64;                         %[deg] Orbital Inclination
orb.raan = 0;                       %[deg] Right Ascension of Ascending Node
orb.e = 0;                          %[] Orbital eccentricity
orb.w = 0;                          %[deg] Argument of Perigee
orb.t0 = 0;                         %[s] Epoch time
orb.a = 6.821e6;                    %[m] semi-major axis
orb.T = 2*pi* sqrt(orb.a^3/Eth.mu);

% Compute Rotation Matrix C_PG (ECI to perifocal)
C_PG = FrameG2P(orb.i, orb.raan, orb.w);
C_GP = C_PG';

%% ---Initial conditions & Controller Setup--- %%
orbit_count = 3;
eps_0 = [0,0,0]';                           % Initial condition
eta_0 = sqrt(1-eps_0'*eps_0);               % Quaternion constraint
theta_dot_0 = [0.01, -0.01, -0.01]'*5;
I = diag([0.0111, 0.0111, 0.0044]);
dt = 0.1;

%% Code Setup
t_final = ceil(orb.T * orbit_count);     % ceiling to cover entire mission;
t_ary = 0:dt:t_final;
len_t = length(t_ary);   
I_inv = inv(I);

%% Attitude Control Actuators
% Magnetic Torquer
MagTorq.n = 1;                % [Turns] turns of each torquer
MagTorq.S = 0.238;   % [m^2] vector area of each torquer
MagTorq.R = 23.15;                % [ohm] Resistance
MagTorq.cur = zeros(3,len_t);   % [A] Current
MagTorq.P = zeros(3,len_t);     % [W] Power

% Magnetic Hysteresis
% ---
% [INSERT HYSTERESIS PROPERTIES HERE]
% ---

%% RK4 Forward Integration
X_qua = zeros(4,len_t);         % Quaternion state vector, [epsilon; eta]
X_ang = zeros(6,len_t);         % [theta; theta_dot], using small angle approximation
omega = zeros(3,len_t);         % Angular velocity         
u_ct_ary = zeros(3,len_t);      % Magnetic dipole of permanent magnet
torq.mag = zeros(3,len_t);      % 3xN array of torque by permanent magnet (in body-axis, no z-axis torque)
torq.hys = zeros(3,len_t);      % 3xN array of hysteresis torque (should be no z-axis torque)
torq.avi = zeros(3,len_t);      % Magnetic dipole from avionics/ steppers
torq.grav = zeros(3,len_t);     % Gravitational torque
torq.sum = zeros(3,len_t);      % Sum of all torques


%% Performance Indices
PerfParam.mag_tau = 0;
PerfParam.omg = 0;
PerfParam.phi = 0;

% Initial Conditions
X_qua_0 = [eps_0;eta_0];
omega(:,1) = theta_dot_0;
X_ang_0 = [2*eps_0 ; theta_dot_0];
X_qua(:,1) = X_qua_0;
X_ang(:,1) = X_ang_0;

% Other parameter setups
E_prev = 0;

for indx1 = 2:len_t
    X_qua_prev = X_qua(:,indx1-1);
    X_ang_prev = X_ang(:,indx1-1);
    omega_prev = omega(:,indx1-1);
    t_prev = t_ary(indx1-1);
        
    % --- Compute Continous Torque --- %
    u_ct = [0 0 0.5];       % MODIFY: magnetic dipole of the permanent magnet
    
    % --- Forward Integrate Attitude --- %    
    [X_qua_now, omega_now, torques, E_new, PerfParam] =...
        RK4_NonLinear(dt,t_prev, X_qua_prev, omega_prev, u_ct, orb,Eth,I,I_inv,C_GP,E_prev,PerfParam);
    
    u_ct_ary(:,indx1) = u_ct;
    torq.mag(:,indx1) = torques.mag;
    torq.avi(:,indx1) = torques.avi;
    torq.grav(:,indx1) = torques.gravgrad;
    torq.sum(:,indx1) = torques.sum;
    
    X_qua(:,indx1) = X_qua_now;
    omega(:,indx1) = omega_now;
    X_ang(1:3,indx1) = 2*X_qua_now(1:3);
    X_ang(4:6,indx1) = omega_now;
    E_prev = E_new;
    
end

%% Compute Performance Parameters
func_eval = @(integral) sqrt(integral/(orbit_count*orb.T));
PerfParam.mag_tau = func_eval(PerfParam.mag_tau);
PerfParam.omg = func_eval(PerfParam.omg);
PerfParam.phi = func_eval(PerfParam.phi);
X_end = X_ang(:,len_t);                         % X(t_f)

fprintf('done\n');

%% Plotter
t_axis = t_ary./orb.T;
% ReportPlotter
Plotter3(11,t_axis,X_ang(1:3,:),'\theta','rad',X_ang(4:6,:),'\omega','rad/s')

figure(6)
subplot(3,1,1)
plot(t_axis,torq.mag(1,:));
title('Permenant magnet torque in 3 axis')
ylabel('X-axis torque (Nm)')
subplot(3,1,2)
plot(t_axis,torq.mag(2,:));
ylabel('Y-axis torque (Nm)')
subplot(3,1,3)
plot(t_axis,torq.mag(3,:));
ylabel('Z-axis torque (Nm)')

Plotter3(4,t_axis,torq.avi,'\tau_a_v_i ','Nm',torq.grav,'\tau_g_r_a_v ','Nm')
Plotter3(1,t_axis,X_qua(1:3,:),'\epsilon','',omega,'\omega','rad/s')

%% Print performance indices to command prompt
fprintf('||tau_mag||_%dT = %.3e\n',orbit_count, PerfParam.mag_tau);
fprintf('|| omega ||_%dT = %.3e\n',orbit_count, PerfParam.omg);
fprintf('||  phi  ||_%dT = %.3e\n',orbit_count, PerfParam.phi);

%% Quaternion checker
eps = X_qua_now(1:3);
eta = X_qua_now(4);
checker = eps'*eps + eta^2;
fprintf('final quaternion check: %.10f\n',checker);

