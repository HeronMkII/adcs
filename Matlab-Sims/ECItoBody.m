function [Rbi] = ECItoBody(q)
%Given the quaternion expression of the body in the ECI reference frame in 
%form (theta;a1;a2;a3, creates the rotation cosine matrix from ECI to Body.

Rbi = [q(2)^2-q(3)^2-q(4)^2+q(1)^2, 2*(q(2)*q(3)+q(1)*q(4)), 2*(q(2)*q(4)-q(1)*q(3));...
       2*(q(2)*q(3)-q(1)*q(4)), -q(2)^2+q(3)^2-q(4)^2+q(1)^2, 2*(q(3)*q(4)+q(1)*q(2));...
       2*(q(2)*q(4)+q(1)*q(3)), 2*(q(3)*q(4)-q(1)*q(2)), -q(2)^2-q(3)^2+q(4)^2+q(1)^2];

end

