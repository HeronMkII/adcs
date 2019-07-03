function [C_BG] = Qua2CBG(eps,eta)
%QUA2CBG Computes the rotation matrix C_BG
%   Using the quaternion parameters, the rotation matrix between the
%   inertial frame to the body-fixed frame is returned

holder1 = (eta^2 - (eps')*eps) * eye(3);
holder2 = 2*eps*(eps');
holder3 = 2*eta*cross_mat(eps);

C_BG = holder1 + holder2 - holder3;

end

