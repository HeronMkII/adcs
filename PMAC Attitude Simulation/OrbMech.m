function [r_G, E] = OrbMech(orb,Eth,C_GP, t, E_prev)
%ORBMECH models the orbital mechanics of the spacecraft
%   Assuming an ideal, non-perturbed orbit, the position of
%   the spacecraft (in the inertial frame) is determined, using the ISS
%   orbital elements

% Solve for Eccentric Anomaly (E)
E_fun = @(E) E-orb.e*sin(E) - sqrt(Eth.mu/orb.a^3) * (t-orb.t0);
E = fzero(E_fun,E_prev);

% Solve for True Anomaly (theta)
holder1 = sqrt((1+orb.e)/(1-orb.e)) * tan(E/2);
theta = atan(holder1)*2;

% Compute SC position vector in the ECI frame (r_G)
r = orb.a*(1-orb.e^2)/(1+orb.e*cos(theta));
r_G = C_GP * [r*cos(theta); r*sin(theta); 0];

end

