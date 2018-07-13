function [B] = dipole_magfield(r)
%Returns magnetic field in ECF coordinates, assuming a dipole model for the
%Earth. Core is considered a giant magnet.

Me=[0;0;-7.94e22];
mu=4e-7*pi;

B = (mu/(4*pi))*(3*r*sum(Me.*r)/norm(r)^5 - Me/norm(r)^3);

end

