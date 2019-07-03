function [Xdot] = OmegaQuad_dot(X_prev, I,I_inv, torq, h_w)
%OMEGAQUAD_DOT computes the time derivative of the attitude states
%   Attitude states refer to [omega, epsilon, eta]'
    
    w = X_prev(1:3,1);      % omega vector
    q = X_prev(4:7,1);      % quaternion vector [eps; eta]
    Mat = 1/2 * [-cross_mat(w), w; -w' 0];  % SDAC II Chapter4 Eqn(9)

    wdot = I_inv*(torq - cross_mat(w)* (I*w+h_w) );     % obtain time rate of change of omega vector
    qdot = Mat * q;                         % obtain time rate of change of quaternion
    
    Xdot = [wdot;qdot];                     % return time derivative of state vector
    
end

