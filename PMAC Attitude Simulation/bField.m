function [b_G] = bField(pos_G)
%BFIELD Models Earth's magnetic field (for the magnetic-torquers)
%   This function admits the cartesian position of the s/c in the ECI frame
%   (G-frame), and returns the B-vector experienced by the s/c
    x = pos_G(1); y = pos_G(2); z = pos_G(3);
    r = norm(pos_G);
    B_0 = -8e15;            % T m^3 (from Prof. Damaren/ online source)
    b_G = B_0/r^5 * [3*x*z,...
                    3*y*z,...
                    2*z^2-x^2-y^2]';
end

