function [X_ang_dot, E_new,b_G_cross] = Ang_diff(A_ct, X_ang_prev, u_ct, t_prev, orb, Eth, I_inv, C_GP,E_prev)
%ANG_DIFF computes the rate of change of the angular state vector
%   Implements continous attitude dynamics, equation (6a) of Behrad's paper

% Control Input
[B_ct,E_new,b_G_cross] = Bct_dynamics(I_inv, orb,Eth,C_GP,t_prev,E_prev);

X_ang_dot = A_ct * X_ang_prev + B_ct * u_ct; 
end

