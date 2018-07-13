function [H] = dipole_magstrength(r)
%Returns magnetic field strength in ECF coordinates, assuming a dipole model for the
%Earth. Core is considered a giant magnet.

Me=[0;0;-7.94e22];

H = (1/(4*pi))*(3*r*Me'*r/(norm(r))^5 - Me/norm(r)^3);

end

