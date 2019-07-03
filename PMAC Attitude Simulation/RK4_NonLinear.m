function [qua_now, omega_now, torques, E_new, PerfParam] = RK4_NonLinear(dt,t_prev, qua_prev, omega_prev, u_ct,orb,Eth,I,I_inv,C_GP,E_prev,PerfParam)
%RK4_NONLINEAR integrates the system of non-linear differential equations
%   Uses 4th order Runge-Kuatta method for the attitude propagation model

    eps_prev = qua_prev(1:3);
    eta_prev = qua_prev(4);
    
%% Compute torque applied
    %% Magnetic Torque
    C_BG_prev = Qua2CBG(eps_prev,eta_prev);
    [r_G_prev, E_new] = OrbMech(orb,Eth,C_GP, t_prev,E_prev);       % Keplerian orbit assumed (function of time only)
    b_G = bField(r_G_prev);
    torques.mag = cross_mat(u_ct) * C_BG_prev * b_G;
    
    %% Gravitational Gradient
    r_B = C_BG_prev * r_G_prev;
    temp1 = 3 * Eth.mu / (norm(r_B))^5;
    temp2 = cross_mat(r_B) * I * r_B;
    torques.gravgrad = temp1 * temp2;

    %% Avionics Disturbance
    m_avi = [0.01; 0.01; 0.01];
    torques.avi = cross_mat(m_avi) * C_BG_prev * b_G;
    
    %% Consolidate Torques
    tau = torques.mag + torques.gravgrad + torques.avi;         % + torques.hys
    torques.sum = tau;
   
    
%% Compute non-linear equations
    h_w = [0;0;0];      % No reaction wheel present in PMAC
    % Step 1
    X_prev = [omega_prev; eps_prev; eta_prev];
    Xdot1 = OmegaQuad_dot(X_prev, I,I_inv, tau, h_w);
    % Step 2
    X_h1 = X_prev + Xdot1*dt/2;
    Xdot2 = OmegaQuad_dot(X_h1, I,I_inv, tau, h_w);
    % Step 3
    X_h2 = X_prev + Xdot2*dt/2;
    Xdot3 = OmegaQuad_dot(X_h2, I,I_inv, tau, h_w);
    % Step 4
    X_f = X_prev + Xdot3*dt;
    Xdot4 = OmegaQuad_dot(X_f, I,I_inv, tau, h_w);
    % Consolidate
    X_now = X_prev + (Xdot1 + 2*Xdot2 + 2*Xdot3 + Xdot4)/6 * dt;
    omega_now = X_now(1:3);
    qua_now = X_now(4:7);
    
%% Performance Parameters
    PerfParam.mag_tau = PerfParam.mag_tau + torques.mag' * torques.mag * dt;
    PerfParam.omg = PerfParam.omg + omega_now' * omega_now * dt;
    phi = acos(0.5*(trace(C_BG_prev)-1));
    PerfParam.phi = PerfParam.phi + phi^2 * dt;
    
end
    
