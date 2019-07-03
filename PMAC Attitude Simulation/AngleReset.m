function [new_angle] = AngleReset(angle)
%ANGLERESET resets the angle value if it is greater than 2pi or smaller
%than 0
%   Re-computes angle value if outside {0,2pi}

if angle >= 2*pi
    new_angle = angle-2*pi;
elseif angle < 0
    new_angle = angle + 2*pi;
else
    new_angle = angle;
end

end

